#ifndef HOOMD_NPH_META_MOVE_CPP
#define HOOMD_NPH_META_MOVE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdNPHMetaMove.h"

#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/simulation/Simulation.h>
#include <mcMd/potentials/pair/McPairPotential.h>
#include <util/ensembles/BoundaryEnsemble.h>
#include <util/crystal/PointGroup.h>
#include <util/crystal/PointSymmetry.h>
#include <util/space/Vector.h>
#include <util/space/Dimension.h>
#include <util/math/Constants.h>
#include <util/math/feq.h>
#include "HoomdNPHMetaMove.cuh"
#include "cudacpu_vector_types.h"
#include "cudacpu_vector_functions.h"

namespace McMd
{

   using namespace Util;

   /*
   * Constructor
   */
   HoomdNPHMetaMove::HoomdNPHMetaMove(McSystem& system) :
      HoomdMove(system), W_(0.0)
   {
   }

   /*
   * Destructor.
   */
   HoomdNPHMetaMove::~HoomdNPHMetaMove()
   {}

   /*
   * Read parameters
   */
   void HoomdNPHMetaMove::readParameters(std::istream& in)
   {
      if ((! system().boundaryEnsemble().isIsobaric()) || (!energyEnsemble().isIsothermal()))
         UTIL_THROW("Must be in isothermal-isobaric ensemble.");

      readProbability(in);
      read<int>(in, "nStep", nStep_);
      read<double>(in, "dt", dt_);

      read<std::string>(in, "mode", modeIn_);
      if (modeIn_ == "cubic")
         integrationMode_ = TwoStepNPHGPU::cubic;
      else if (modeIn_ == "orthorhombic")
         integrationMode_ = TwoStepNPHGPU::orthorhombic;
      else if (modeIn_ == "tetragonal")
         integrationMode_ = TwoStepNPHGPU::tetragonal;
      else
         UTIL_THROW("Unsupported integration mode.");

      read<double>(in,"W", W_);
      read<double>(in, "skin", skin_);
      char* env;
      if ((env = getenv("OMPI_COMM_WORLD_LOCAL_RANK")) != NULL) {
         GPUId_ = atoi(env);
      } else {
         GPUId_ = -1;
      }

      read<std::string>(in, "phase", phase_);
      if (phase_ != "disordered" && phase_ != "ordered")
         UTIL_THROW("Unsupported phase.");

      read<double>(in, "maxSF", maxSF_);

      nAtomType_ = system().simulation().nAtomType();

      read<int>(in, "nMode", nMode_);
      modes_.allocate(nMode_, nAtomType_);
      readDMatrix<double>(in, "modes", modes_, nMode_, nAtomType_);
      read<int>(in, "hMax", hMax_);
      read<Util::LatticeSystem>(in, "lattice", lattice_);

      // Allocate wavevectors arrays
      nWave_     = (2*hMax_ + 1)*(2*hMax_ + 1)*(2*hMax_ + 1);
      waveIntVectors_.allocate(nWave_);
      waveVectors_.allocate(nWave_);
      fourierModes_.allocate(nWave_, nMode_);
      structureFactors_.allocate(nWave_, nMode_);

      int i, j, h, k, l, m;
      IntVector g;

      // Cubic Symmetry
      if (lattice_ == Cubic) {

         nStar_ = (hMax_ + 1 )*(hMax_ + 2)*(hMax_ + 3)/6;
         starIds_.allocate(nStar_);
         starSizes_.allocate(nStar_);

         // Create cubic point group
         PointGroup group;
         PointSymmetry a, b, c;

         a.R(0,1) =  1;
         a.R(1,0) =  1;
         a.R(2,2) =  1;

         b.R(0,0) = -1;
         b.R(1,1) =  1;
         b.R(2,2) =  1;

         c.R(0,1) =  1;
         c.R(1,2) =  1;
         c.R(2,0) =  1;

         group.add(c);
         group.add(b);
         group.add(a);
         group.makeCompleteGroup();

         // Create grid of wavevectors
         FSArray<IntVector, 48> star;
         i = 0;
         j = 0;
         for (h = 0; h <= hMax_; ++h) {
            g[0] = h;
            for (k = 0; k <= h; ++k) {
               g[1] = k;
               for (l = 0; l <= k; ++l) {
                  g[2] = l;
                  starIds_[i] = j;
                  group.makeStar(g, star);
                  starSizes_[i] = star.size();
                  for (m = 0; m < star.size(); ++m) {
                     waveIntVectors_[j] = star[m];
                     ++j;
                  }
                  ++i;
               }
            }
         }
         if (i != nStar_) {
            UTIL_THROW("Error");
         }
         if (j != nWave_) {
            UTIL_THROW("Error");
         }
      } else if (lattice_ == Tetragonal) {

         nStar_ = (hMax_ + 1 )*(hMax_ + 1)*(hMax_ + 2)/2;
         starIds_.allocate(nStar_);
         starSizes_.allocate(nStar_);
         
         // Create tetragonal point group
         PointGroup group;
         PointSymmetry a, b, c;

         a.R(0,0) =  1;
         a.R(1,2) =  1;
         a.R(2,1) =  1;

         b.R(0,0) =  -1;
         b.R(1,1) =  1;
         b.R(2,2) =  1;

         c.R(0,0) =  1;
         c.R(1,1) =  -1;
         c.R(2,2) =  1;

         group.add(c);
         group.add(b);
         group.add(a);
         group.makeCompleteGroup();
         
         // Create grid of wavevectors
         FSArray<IntVector, 16> star;
         i = 0;
         j = 0;
         for (h = 0; h <= hMax_; ++h) {
            g[0] = h;
            for (k = 0; k <= hMax_; ++k) {
               g[1] = k;
               for (l = 0; l <= k; ++l) {
                  g[2] = l;
                  starIds_[i] = j;
                  group.makeStar(g, star);
                  starSizes_[i] = star.size();
                  for (m = 0; m < star.size(); ++m) {
                     waveIntVectors_[j] = star[m];
                     ++j;
                  }
                  ++i;
               }
            }
         }
         if (i != nStar_) {
            UTIL_THROW("Error");
         }
         if (j != nWave_) {
            UTIL_THROW("Error");
         }
      }
      
      // Clear accumulators
      for (i = 0; i < nWave_; ++i) {
         for (j = 0; j < nMode_; ++j) {
            structureFactors_(i, j) = 0.0;
         }
      }

      toImposeConstrain_ = false;
      setConstrain_ = false;

      //read<int>(in, "GPUId", GPUId_);
      // create HOOMD execution configuration
      executionConfigurationSPtr_ = 
        boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::GPU,GPUId_));

      // set CUDA error checking
      #ifdef UTIL_DEBUG
      executionConfigurationSPtr_->setCUDAErrorChecking(true);
      #endif

      // subscribe to molecule set notifications
      system().subscribeMoleculeSetChange(*this);

      #if defined(HOOMD_DEVEL) && defined(MCMD_LINK)
      // subscribe to link notifications
      system().linkMaster().Notifier<LinkAddEvent>::registerObserver(*this);
      system().linkMaster().Notifier<LinkResetEvent>::registerObserver(*this);
      system().linkMaster().Notifier<LinkRemoveEvent>::registerObserver(*this);
      #endif 

   }
 
   void HoomdNPHMetaMove::createIntegrator()     
   {

      // create NVE Integrator 
      integratorSPtr_ = boost::shared_ptr<IntegratorTwoStep>(new IntegratorTwoStep(systemDefinitionSPtr_,dt_));

      boost::shared_ptr<Variant> variantPSPtr(new VariantConst(system().boundaryEnsemble().pressure()));
      // create IntegrationMethod
      twoStepNPHGPUSPtr_ = boost::shared_ptr<TwoStepNPHGPU>(new TwoStepNPHGPU(systemDefinitionSPtr_,groupAllSPtr_, thermoSPtr_, W_, variantPSPtr, integrationMode_,""));

      integratorSPtr_->addIntegrationMethod(twoStepNPHGPUSPtr_);

      // register pair and bond forces in Integrator
      integratorSPtr_->addForceCompute(pairForceSPtr_);
      integratorSPtr_->addForceCompute(bondForceSPtr_);
      #ifdef INTER_EXTERNAL
      if (implementExternalPotential_)
         integratorSPtr_->addForceCompute(externalForceSPtr_);
      #endif


      // set flags for thermodynamic quantities requested by integrator
      // additionally, we require computation of potential energy
      PDataFlags peFlag;
      peFlag[pdata_flag::potential_energy] = 1;
      particleDataSPtr_->setFlags(integratorSPtr_->getRequestedPDataFlags() | peFlag);
   }

   /*
   * Generate, attempt and accept or reject a Hybrid MD/MC move.
   */
   bool HoomdNPHMetaMove::move()
   {
      if ((!HoomdIsInitialized_) || moleculeSetHasChanged_) {
         initSimulation();
         moleculeSetHasChanged_ = false;
      }

      if ( !setConstrain_ && toImposeConstrain_) {
         Boundary initBoundary = system().boundary();
         constrainLengths_ = initBoundary.lengths();
         setConstrain_ = true;
      }

      // We need to create the Integrator every time since we are starting
      // with new coordinates, velocities etc.
      // this does not seem to incur a significant performance decrease
      createIntegrator();

      System::MoleculeIterator molIter;
      Molecule::AtomIterator   atomIter;
      int nSpec = simulation().nSpecies();

      // Increment counter for attempted moves
      incrementNAttempt();

      // set hoomd simulation box
      BoxDim box;
      Boundary boundary = system().boundary();
      lengths_ = boundary.lengths();
      const Scalar3 boundaryLengths = make_scalar3(lengths_[0], lengths_[1], lengths_[2]);
      box.setL(boundaryLengths);

      particleDataSPtr_->setGlobalBoxL(boundaryLengths);
      {
      // copy atom coordinates into hoomd
      ArrayHandle<Scalar4> h_pos(particleDataSPtr_->getPositions(), access_location::host, access_mode::readwrite);
      ArrayHandle<Scalar4> h_vel(particleDataSPtr_->getVelocities(), access_location::host, access_mode::readwrite);
      ArrayHandle<unsigned int> h_tag(particleDataSPtr_->getTags(), access_location::host, access_mode::readwrite);
      ArrayHandle<unsigned int> h_rtag(particleDataSPtr_->getRTags(), access_location::host, access_mode::readwrite);
      int nind = 0;
      for (int iSpec =0; iSpec < nSpec; ++iSpec) {
         system().begin(iSpec, molIter);
         for ( ; molIter.notEnd(); ++ molIter) {
            for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
               nind = nind + 1;
               unsigned int idx = (unsigned int) atomIter->id();
               Vector& pos = atomIter->position();
               h_pos.data[idx].x = pos[0] - lengths_[0]/2.;
               h_pos.data[idx].y = pos[1] - lengths_[1]/2.;
               h_pos.data[idx].z = pos[2] - lengths_[2]/2.;
    
               int type = atomIter->typeId();
               h_vel.data[idx].w = simulation().atomType(type).mass();
               h_pos.data[idx].w = __int_as_scalar(type);
               h_tag.data[idx] = idx;
               h_rtag.data[idx] = idx; 
               /*
               if(nind == 1) {
                  std::cout << " index is " << atomIter->id() << std::endl;
                  std::cout << "position is " << pos << std::endl;
                  std::cout << "lengths is " << lengths_ << std::endl;
                  std::cout << "hpos is " << h_pos.data[idx].x << "  " << h_pos.data[idx].y << " " << h_pos.data[idx].z << std::endl;
               }*/
            }
         }
      }

      // Generate random velocities
      generateRandomVelocities(h_vel);
      }

      // Done with the data

      // generate integrator variables from a Gaussian distribution
      double etax = 0.0;
      double etay = 0.0;
      double etaz = 0.0;
      Random& random = simulation().random();
      double temp = energyEnsemble().temperature();
      if (integrationMode_ == TwoStepNPHGPU::cubic) {
         // one degree of freedom
         // barostat_energy = 1/2 (1/W) eta_x^2
         double sigma = sqrt(temp/W_);
         etax = sigma*random.gaussian();
      } else if (integrationMode_ == TwoStepNPHGPU::tetragonal) {
         // two degrees of freedom
         // barostat_energy = 1/2 (1/W) eta_x^2 + 1/2 (1/(2W)) eta_y^2
         double sigma1 = sqrt(temp/W_);
         etax = sigma1*random.gaussian();
         double sigma2 = sqrt(temp/W_/2.0);
         etay = sigma2*random.gaussian();
      } else if (integrationMode_ == TwoStepNPHGPU::orthorhombic) {
         // three degrees of freedom 
         // barostat_energy = 1/2 (1/W) (eta_x^2 + eta_y^2 + eta_z^2)
         double sigma = sqrt(temp/W_);
         etax = sigma*random.gaussian();
         etay = sigma*random.gaussian();
         etaz = sigma*random.gaussian();
      } 
      twoStepNPHGPUSPtr_->setEta(etax,etay,etaz);
      
      // Notify that the particle order might have changed
      particleDataSPtr_->notifyParticleSort();

      // Calculate oldH (the conserved quantity of the Andersen barostat)
      double volume = lengths_[0]*lengths_[1]*lengths_[2];

      // Initialize integrator (calculate forces and potential energy for step 0)
      integratorSPtr_->prepRun(0);

      // H = U + PV + barostat_energy
      thermoSPtr_->compute(0);
      double oldH = thermoSPtr_->getLogValue("kinetic_energy",0);
      oldH += thermoSPtr_->getLogValue("potential_energy",0);
      oldH += system().boundaryEnsemble().pressure() * volume;
      oldH += integratorSPtr_->getLogValue("nph_barostat_energy",0);

      // Integrate nStep_ steps forward
      for (int iStep = 0; iStep < nStep_; ++iStep) {
         integratorSPtr_->update(iStep);
         //std::cout << "energy at iStep is " << externalForceSPtr_->getLogValue("external_periodic_energy", iStep) << std::endl; 

         // do we need to sort the particles?
         // do not sort at time step 0 to speed up short runs
         if (! (iStep % sorterPeriod_) && iStep)
           sorterSPtr_->update(iStep);
      }

      Vector newLengths;
      box = particleDataSPtr_->getBox();
      newLengths[0] = (box.getL()).x;
      newLengths[1] = (box.getL()).y;
      newLengths[2] = (box.getL()).z;
      volume = newLengths[0]*newLengths[1]*newLengths[2];
 
      // Calculate new value of the conserved quantity
      thermoSPtr_->compute(nStep_);
      double newH = thermoSPtr_->getLogValue("kinetic_energy",nStep_);
      newH += thermoSPtr_->getLogValue("potential_energy",nStep_);
      newH += system().boundaryEnsemble().pressure() * volume;
      newH += integratorSPtr_->getLogValue("nph_barostat_energy",nStep_);

      ArrayHandle<Scalar4> h_pos_integration(particleDataSPtr_->getPositions(), access_location::host, access_mode::read);
      ArrayHandle<unsigned int> h_tag_integration(particleDataSPtr_->getTags(), access_location::host, access_mode::read);
      ArrayHandle<unsigned int> h_rtag_integration(particleDataSPtr_->getRTags(), access_location::host, access_mode::read);

      Vector  position;
      int  typeId;

      // Calculate wavevectors
      for (int i = 0; i < nWave_; ++i) {
         waveVectors_[i] = Vector::Zero;
         for (int j = 0; j < Dimension; ++j) {
            waveVectors_[i][j] = 2.0*Constants::Pi*waveIntVectors_[i][j]/newLengths[j];
         }
      }

      // allocate temporary arrays
      float3 *h_wave_vectors = new float3[nWave_];
      float3 *h_position = new float3[system().simulation().atomCapacity()];
      int *h_type = new int[system().simulation().atomCapacity()];
      float *h_mode = new float[nMode_*nAtomType_];
      float *h_sq = new float[nMode_*nWave_];

      // Load wave vectors
      for (int i = 0; i < nWave_; i++)
         h_wave_vectors[i] = make_float3((float) waveVectors_[i][0],
                                         (float) waveVectors_[i][1],
                                         (float) waveVectors_[i][2]);

      // Load atoms into temporary storage
      int nAtom = 0;
      for (int iSpec = 0; iSpec < nSpec; ++iSpec) {
         system().begin(iSpec, molIter);
         for ( ; molIter.notEnd(); ++molIter) {
            molIter->begin(atomIter);
            for ( ; atomIter.notEnd(); ++atomIter) {
               unsigned int idx = h_rtag_integration.data[atomIter->id()];
               typeId   = atomIter->typeId();
               position  = Vector(h_pos_integration.data[idx].x+newLengths[0]/2.,
                                                h_pos_integration.data[idx].y+newLengths[1]/2.,
                                                h_pos_integration.data[idx].z+newLengths[2]/2.);
               h_position[nAtom] = make_float3(position[0],position[1],position[2]);
               h_type[nAtom] = typeId;
               nAtom++;
            }

         }
      }

      // Load modes
      for (int i = 0; i < nMode_; ++i) {
         for (int j = 0; j < nAtomType_; ++j) {
               h_mode[nAtomType_*i + j] = modes_(i, j);
         }
      }

      // Calculate structure factors
      int res = gpu_sample_structure_factor(nWave_,
                                     h_wave_vectors,
                                     nAtom,
                                     h_position,
                                     h_type,
                                     nAtomType_,
                                     nMode_,
                                     h_mode,
                                     h_sq,
                                     volume
                                     );

      if (res)
          UTIL_THROW("Error updating structure factor.");

      bool accept = true;

      double maxValue;

      // increment structure factors
      for (int j = 0; j < nMode_; ++j) {
         maxValue = 0.0;
         for (int i = 0; i < nWave_; ++i) {
            if ((double) h_sq[j*nWave_+i] >= maxValue) {
               maxValue = (double) h_sq[j*nWave_+i];
            }
         }
      }

      if (phase_ == "disordered") {
         if (maxValue > maxSF_)
            accept = false;
      } else if (phase_ == "ordered") {
         if (maxValue < maxSF_)
            accept = false;
      } 

      // Decide whether to accept or reject
      if (accept != false) {
      if (integrationMode_ == TwoStepNPHGPU::tetragonal) {
         if (!setConstrain_) {
            double  newAspectRatio, aspectRatioParam;
            newAspectRatio = double(newLengths[0]/newLengths[1]);
            if ( newAspectRatio > 1.4 || newAspectRatio < 0.8) {
               accept = false;
            } else {
               accept = random.metropolis( boltzmann(newH-oldH) );
            }
         } else
         if (setConstrain_) {
            double diffLx, diffLy;
            if ( newLengths[0] >= constrainLengths_[0] ) {
               diffLx = newLengths[0] - constrainLengths_[0];
            } else {
               diffLx = -(newLengths[0] - constrainLengths_[0]);
            }
            if ( newLengths[1] >= constrainLengths_[1] ) {
               diffLy = newLengths[1] - constrainLengths_[1];
            } else {
               diffLy = -(newLengths[1] - constrainLengths_[1]);
            }
            if ( diffLx > 0.5 || diffLy > 0.5 ) {
               accept = false;
            } else {
               accept = random.metropolis( boltzmann(newH-oldH) );
            }
         }
      } else {
            accept = random.metropolis( boltzmann(newH-oldH) );
      }
      }

      if (accept) {
         // read back new boundary
         box = particleDataSPtr_->getBox();
         lengths_ = newLengths;
         system().boundary().setOrthorhombic(lengths_);
         
         // read back integrated positions
         ArrayHandle<Scalar4> h_pos(particleDataSPtr_->getPositions(), access_location::host, access_mode::read);
         ArrayHandle<Scalar4> h_vel(particleDataSPtr_->getVelocities(), access_location::host, access_mode::read);
         ArrayHandle<unsigned int> h_tag(particleDataSPtr_->getTags(), access_location::host, access_mode::read);
         ArrayHandle<unsigned int> h_rtag(particleDataSPtr_->getRTags(), access_location::host, access_mode::read);
         //int nind = 0; 
         for (int iSpec = 0; iSpec < nSpec; ++iSpec) {
            system().begin(iSpec, molIter);
            for ( ; molIter.notEnd(); ++molIter) {
               for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
                  //nind = nind + 1;
                  unsigned int idx = h_rtag.data[atomIter->id()]; 
                  atomIter->position() = Vector(h_pos.data[idx].x+lengths_[0]/2.,
                                                h_pos.data[idx].y+lengths_[1]/2.,
                                                h_pos.data[idx].z+lengths_[2]/2.);
                  //if(nind == 85408) {
                    // std::cout << "accept hpos is " << h_pos.data[idx].x << std::endl;
                    // std::cout << "accept position is " << atomIter->position() << std::endl;
                  //}
               }
            }
         }
         // Done with data
         system().pairPotential().buildCellList();
         incrementNAccept();
      } else {
          // not accepted, do nothing
         // read back integrated positions
         int nind = 0;
         for (int iSpec = 0; iSpec < nSpec; ++iSpec) {
            system().begin(iSpec, molIter);
            for ( ; molIter.notEnd(); ++molIter) {
               for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
                  nind = nind + 1;
                  /*
                  if(nind == 1) {
                     std::cout << "reject index is " << atomIter->id() << std::endl;
                     std::cout << "reject position is " << atomIter->position() << std::endl;
                     std::cout << "lengths is " << lengths_ << std::endl;
                  }*/
               }
            }
         }

      }

//      std::cout << "Lx = " << lengths_[0] << " Ly = " << lengths_[1] << " Lz = " << lengths_[2] << std::endl;
      return accept;
   }
 
}
#endif
