#ifdef INTER_EXTERNAL
#ifndef HOOMD_ORDERING_EXTERNAL_H
#define HOOMD_ORDERING_EXTERNAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdExternal.h"

#include <hoomd/PotentialExternal.h>
#include <hoomd/EvaluatorExternalPeriodic.h>
#include <hoomd/AllDriverPotentialExternalGPU.cuh>

namespace McMd
{
   extern char classNameHoomdOrdering[];

   /**
   * A potential encapsulating the HOOMD Periodic evaluator
   *
   * \ingroup External_Module
   */
   class HoomdOrderingExternal : public HoomdExternal< EvaluatorExternalOrdering,
      gpu_compute_ordering_forces, classNameHoomdOrdering >
   {
   
   public:

      /**
      * Default constructor.
      */
      HoomdOrderingExternal();

      /**
      * Copy constructor
      */
      HoomdOrderingExternal(const HoomdOrderingExternal& other);

      /**
      * read parameters from file
      *
      * \param in input stream
      */
      void readParameters(std::istream &in);

      /**
      * Set external potential parameter
      *
      */
      void setExternalParameter(double externalParameter);

      /**
      * returns external potential parameter
      */
      double externalParameter() const;

      /**
      * returns the class name
      */
      std::string className() const;

   private:

      /// per-type prefactor of potential
      DArray<double> prefactor_;

      /// external parameter
      double externalParameter_;

      /// Number of reciprocal lattice vectors
      int  nWaveVectors_;

      /// Array of Miller index IntVectors for the reciprocal lattice vectors.
      DArray<IntVector>  waveIntVectors_;

      /// Array of floating point reciprocal lattice vectors.
      DArray<Vector>  waveVectors_;

      double interfaceWidth_;

      int periodicity_;
   };
  
}

#endif
#endif
