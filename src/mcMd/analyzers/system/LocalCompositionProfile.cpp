#ifndef INTER_NOPAIR
#ifndef MCMD_LOCAL_COMPOSITION_PROFILE_CPP
#define MCMD_LOCAL_COMPOSITION_PROFILE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "LocalCompositionProfile.h"        // class header

#include <mcMd/simulation/Simulation.h>
#include <mcMd/simulation/McMd_mpi.h>
#include <util/boundary/Boundary.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <util/containers/DMatrix.h>
#include <util/space/Dimension.h>
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <util/format/Dbl.h>
#include <util/archives/Serializable_includes.h>

#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   LocalCompositionProfile::LocalCompositionProfile(System& system)
    : SystemAnalyzer<System>(system),
      outputFile_(),
      isFirstStep_(true),
      isInitialized_(false)
   {
      setClassName("LocalCompositionProfile"); 
   }


   /*
   * Read parameters and initialize.
   */
   void LocalCompositionProfile::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      read<double>(in, "rCutOff", rCutOff_);
      read<int>(in, "nBins", nBins_);
  
      nAtomType_ = system().simulation().nAtomType();
      accumulators_.allocate(nAtomType_, nAtomType_);
      currentAccumulators_.allocate(nAtomType_, nAtomType_);
      logFiles_.allocate(nAtomType_, nAtomType_);
      
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void LocalCompositionProfile::loadParameters(Serializable::IArchive& ar)
   {
      Analyzer::loadParameters(ar);  
      ar & accumulators_;
      ar & rCutOff_;
      ar & nSample_;
      ar & nAtomType_;
      ar & nBins_;
      ar & isFirstStep_;

      if (nAtomType_ != system().simulation().nAtomType()) {
         UTIL_THROW("Inconsistent values for nAtomType_");
      }

      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void LocalCompositionProfile::save(Serializable::OArchive& ar)
   {  ar & *this; }

   /*
   * Clear accumulator.
   */
   void LocalCompositionProfile::setup() 
   {  
      double min, max;
      min = 0.0;
      max = 1.0;

      for (int i=0; i < nAtomType_; ++i) {
        for ( int j=0; j < nAtomType_; ++j) {
          accumulators_(i, j).setParam(min, max, nBins_);
          currentAccumulators_(i, j).setParam(min, max, nBins_);
        }
      }

      // Clear accumulators
      for (int i = 0; i < nAtomType_; ++i){
        for (int j = 0; j < nAtomType_; ++j) {
          accumulators_(i, j).clear();
        }
      }
      nSample_ = 0;
   }

 
   /*
   * Evaluate energy per particle, and add to ensemble. 
   */
   void LocalCompositionProfile::sample(long iStep) 
   {

      if (isAtInterval(iStep)) {;
       
      System::ConstMoleculeIterator molIter1, molIter2;
      Molecule::ConstAtomIterator   atomIter1, atomIter2;
      Vector    r1, r2;
      double    dRsq, dR;
      Boundary* boundaryPtr;
      int       iSpecies1, iSpecies2, nSpecies, m, n;
      DMatrix<double> Ncount;
      Ncount.allocate(nAtomType_, nAtomType_);

      // Clear accumulators for the current timestep
      for (m = 0; m < nAtomType_; ++m){
        for (n = 0; n < nAtomType_; ++n) {
          currentAccumulators_(m, n).clear();
          std::ostringstream oss;
          oss << outputFileName();

          oss << "_type" << m << "_type" << n << ".log";

          fileMaster().openOutputFile(oss.str(), logFiles_(m, n), !isFirstStep_);
        }
      }

      boundaryPtr = &system().boundary();
      nSpecies    = system().simulation().nSpecies();

      // Loop over atom 1
      for (iSpecies1 = 0; iSpecies1 < nSpecies; ++iSpecies1) {
        system().begin(iSpecies1, molIter1);
            
        for ( ; molIter1.notEnd(); ++molIter1) {
          molIter1->begin(atomIter1);
               
          for ( ; atomIter1.notEnd(); ++atomIter1) {
            r1 = atomIter1->position();
            for (m = 0; m < nAtomType_; ++m) {
              for (n = 0; n < nAtomType_; ++n) {
                Ncount(m, n) = 0.0;
              }
            }
            
            for (iSpecies2 = 0; iSpecies2 < nSpecies; ++iSpecies2) {
               system().begin(iSpecies2, molIter2);
               
               for ( ; molIter2.notEnd(); ++molIter2) {
                  if ( &(*molIter2) != &(*molIter1)) {
                     molIter2->begin(atomIter2);
                           
                     for ( ; atomIter2.notEnd(); ++atomIter2) {
                       r2 = atomIter2->position();
                       dRsq = boundaryPtr->distanceSq(r1, r2);
                       dR   = sqrt(dRsq);
                       Ncount(atomIter1->typeId(), atomIter2->typeId()) 
                            += 1/(1 + exp(12*(dR - rCutOff_)));
                     }
                  }
               }
            }
            double num1 = Ncount(atomIter1->typeId(), atomIter1->typeId());
            double denom1 = Ncount(atomIter1->typeId(), atomIter1->typeId()) 
                          + Ncount(atomIter1->typeId(), 1-(atomIter1->typeId()));
            double num2 = Ncount(atomIter1->typeId(), 1-(atomIter1->typeId()));
            double denom2 = Ncount(atomIter1->typeId(), atomIter1->typeId()) 
                          + Ncount(atomIter1->typeId(), 1-(atomIter1->typeId()));
            double fraction1, fraction2;
            if (denom1 == 0.0) {
              fraction1 = 0.0;
            } else {
              fraction1 = num1/denom1;
            }
            if (denom2 == 0.0) {
              fraction2 = 0.0;
            } else {
              fraction2 = num2/denom2;
            }
            currentAccumulators_(atomIter1->typeId(), atomIter1->typeId()).sample(fraction1);
            accumulators_(atomIter1->typeId(), atomIter1->typeId()).sample(fraction1);
            currentAccumulators_(atomIter1->typeId(), 1-(atomIter1->typeId())).sample(fraction2);
            accumulators_(atomIter1->typeId(), 1-(atomIter1->typeId())).sample(fraction2);
          }
        }
      }
      ++nSample_;
      
      for (m = 0; m < nAtomType_; ++m){
        for (n = 0; n < nAtomType_; ++n){
          currentAccumulators_(m, n).output(logFiles_(m, n));
          logFiles_(m, n) << std::endl;
        }
      }
      for (m = 0; m < nAtomType_; ++m){
        for (n = 0; n < nAtomType_; ++n) {
          logFiles_(m, n).close();
        }
      }
      isFirstStep_ = false;
      }
   }

   /*
   * Output results to file after simulation is completed.
   */
   void LocalCompositionProfile::output() 
   { 
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      for (int i = 0; i < nAtomType_; ++i) {
        for (int j = 0; j < nAtomType_; ++j) {
          accumulators_(i, j).output(outputFile_);
          outputFile_ << std::endl;
        }
      }
      outputFile_.close();
   }
   
}
#endif
#endif 
