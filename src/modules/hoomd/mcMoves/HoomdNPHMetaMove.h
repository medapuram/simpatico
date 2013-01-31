#ifndef HOOMD_NPH_META_MOVE_H
#define HOOMD_NPH_META_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdMove.h"

#include <util/crystal/LatticeSystem.h>
#include <hoomd/TwoStepNPHGPU.h>

#include <cmath>
namespace McMd
{

   using namespace Util;

   /**
   * HoomdNPHMetaMove is a hybrid Molecular Dynamics MC move using Hoomd-blue
   * 
   * It samples the NPT (isobaric-isothermal) ensemble by integrating the NPH (isoenthalpic-isobaric,
   * Andersen barostat) equations of motion. The box dimensions are updated after an accepted HoomdNPHMetaMove.
   *
   * \ingroup McMove_Module
   */
   class HoomdNPHMetaMove : public HoomdMove
   {

   public:

      /**
      * Constructor.
      *
      * Constructs a component MdSystem object.
      */
      HoomdNPHMetaMove(McSystem& system);

      /**
      * Destructor.
      */
      ~HoomdNPHMetaMove();

      /**
      * Read nStep, dt, skin, maxNPair from file.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Generate, attempt and accept or reject a move.
      */
      bool move();

      /**
      * Update wavevectors.
      */
      void makeWaveVectors();

   protected:
      /// Create NPH integrator
      void createIntegrator();

      /// Andersen barostat mass
      double W_;
      
      /// anisotropic integration mode (can be cubic, orthorhombic or tetragonal)
      TwoStepNPHGPU::integrationMode integrationMode_; 

      /// integration method
      boost::shared_ptr<TwoStepNPHGPU> twoStepNPHGPUSPtr_;

      /// phase to be maintained
      std::string phase_;
     
      /// limiting value of maximum Structure factor  
      double maxSF_;
      
      /// Number of mode vectors
      int  nMode_;

      /**
      * Array of mode vectors 
      *
      * First index is mode, second is atomType.
      */
      DMatrix<double>  modes_;

      /// Maximum Miller index of wavevectors in grid.
      int   hMax_;

      /// Lattice system used to create stars.
      LatticeSystem   lattice_;

      /**
      * Array of vector of maximum structure factor values. 
      */
      //DArray< std::vector<double> > maximumValue_;

      /**
      * Array of vector of Miller index IntVector with maximum S(q).
      */
      //DArray< std::vector<IntVector> > maximumWaveIntVector_;

      /**
      * Array of vector of magnitudes of waveVector with maximum S(q).
      */
      //DArray< std::vector<double> > maximumQ_;

      /// Number of samples thus far.
      //int  nSample_;

   private:
      
      /// Geometry of simulation cell
      std::string modeIn_;

      /// Number of atom types, copied from Simulation::nAtomType().
      int  nAtomType_;

      /// Number of wavevectors.
      int  nWave_;

      /**
      * Array of Miller index IntVectors for wavevectors.
      */
      DArray<IntVector>  waveIntVectors_;

      /**
      * Array of floating point wave vectors.
      */
      DArray<Vector>  waveVectors_;

      /**
      * Structure factor accumulators. 
      * 
      * First index is wavevector, second index is mode index.
      */
      DMatrix<double> structureFactors_;

      /**
      * Fourier modes of concentration.
      *
      * First index is wavevector, second is atom type.
      */
      DMatrix< std::complex<double> >  fourierModes_;

      /// Array of ids for first wavevector in each star.
      DArray<int>  starIds_;

      /// Array of star sizes.
      DArray<int>  starSizes_;

      /// Number of stars of symmetry related wavevectors.
      int   nStar_;

      bool toImposeConstrain_;
      bool setConstrain_;
      Vector constrainLengths_;
   };

}
#endif
