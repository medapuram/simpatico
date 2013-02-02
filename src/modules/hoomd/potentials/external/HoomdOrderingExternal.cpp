#ifdef INTER_EXTERNAL
#ifndef HOOMD_ORDERING_EXTERNAL_CPP
#define HOOMD_ORDERING_EXTERNAL_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdOrderingExternal.h"


namespace McMd
{
   char classNameHoomdPeriodic[] = "HoomdPeriodicExternal";

   /**
   * Default constructor.
   */
   HoomdOrderingExternal::HoomdOrderingExternal()   
    : HoomdExternal< EvaluatorExternalOrdering, gpu_compute_ordering_forces, classNameHoomdOrdering >()
   {
   }

   /**
   * Copy constructor
   */
   HoomdOrderingExternal::HoomdOrderingExternal(const HoomdOrderingExternal& other)
    : HoomdExternal< EvaluatorExternalOrdering, gpu_compute_ordering_forces, 
          classNameHoomdOrdering >(other)
   {
      externalParameter_ = other.externalParameter_;
      interfaceWidth_ = other.interfaceWidth_;
      periodicity_ = other.periodicity_;
   }

   /**
   * read parameters from file
   */
   void HoomdOrderingExternal::readParameters(std::istream &in)
   {
      // Read parameters
      prefactor_.allocate(nAtomType_);
      readDArray<double>(in, "prefactor", prefactor_, nAtomType_);

      read<double>(in, "externalParameter", externalParameter_);

      read<int>(in, "nWaveVectors", nWaveVectors_);
      waveIntVectors_.allocate(nWaveVectors_);
      readDArray<IntVector>(in, "waveIntVectors", waveIntVectors_, nWaveVectors_);

      waveVectors_.allocate(nWaveVectors_);
      for (int i = 0; i < nWaveVectors_; ++i) {
         waveVectors_[i] = Vector::Zero;
         for (int j = 0; j < Dimension; ++j) {
            Vector dWave  = Vector::Zero;
            dWave[j] = 2.0*M_PI*waveIntVectors_[i][j];
            waveVectors_[i] += dWave;
         }
      }
      read<double>(in, "interfaceWidth", interfaceWidth_);
      read<int>(in, "periodicity", periodicity_);
   
      for (int i = 0; i < nAtomType_; ++i) {
         params_[i].x = __int_as_scalar(perpDirection_);
         params_[i].y = Scalar(prefactor_[i]*externalParameter_);
         params_[i].z = Scalar(width_);
         params_[i].w = __int_as_scalar(periodicity_);
      } 
   }

   /*
   * set external potential parameter
   */
   void HoomdOrderingExternal::setExternalParameter(double externalParameter)
   {
      externalParameter_ = externalParameter;
      for (int i = 0; i < nAtomType_; ++i) {
         params_[i].y = prefactor_[i]*externalParameter_;
      }
   }

   /* 
   * Get external potential interaction strength.
   */
   double HoomdOrderingExternal::externalParameter() const
   {
      return externalParameter_;
   }

   /*
   * return the class name
   */
   std::string HoomdOrderingExternal::className() const
   {
      return "HoomdOrderingExternal";
   }

}

#endif
#endif
