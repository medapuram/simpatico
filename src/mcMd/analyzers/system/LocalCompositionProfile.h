#ifndef INTER_NOPAIR
#ifndef MCMD_LOCAL_COMPOSITION_PROFILE_H
#define MCMD_LOCAL_COMPOSITION_PROFILE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>   // base class template
#include <mcMd/simulation/System.h>          // base template parameter
#include <util/accumulators/Distribution.h>           // member

#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /**
   * LocalCompositionProfile averages of total potential energy.
   *
   * \ingroup McMd_Analyzer_Module
   */
   class LocalCompositionProfile : public SystemAnalyzer<System>
   {
   
   public:

      /**   
      * Constructor.
      *
      * \param system parent McSystem
      */
      LocalCompositionProfile(System& system);

      /**
      * Read parameters and initialize.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to an archive.
      *
      * \param ar saving (output) archive
      */
      virtual void save(Serializable::OArchive& ar);
  
      /**
      * Serialize to/from an archive. 
      * 
      * \param ar      archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /** 
      * Clear accumulator.
      */
      virtual void setup();
   
      /* 
      * Evaluate energy per particle, and add to ensemble. 
      */
      virtual void sample(long iStep);
   
      /**
      * Output results at end of simulation.
      */
      virtual void output();

   private:

      /// Output file stream
      std::ofstream outputFile_;

      double rCutOff_;

      /// Average object - statistical accumulator
      DMatrix<Distribution>  accumulators_;

      DMatrix<Distribution>  currentAccumulators_;

      /// Number of samples thus far.
      int  nSample_;
      
      /// Number of atom types, copied from Simulation::nAtomType().
      int  nAtomType_;
      
      /// Number of bins for density profile
      int nBins_;
      
      /// True if this is the first step
      bool isFirstStep_;

      DMatrix<std::ofstream> logFiles_;

      // Has readParam been called?
      bool isInitialized_;

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void LocalCompositionProfile::serialize(Archive& ar, const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & accumulators_;
   }

}
#endif 
#endif 
