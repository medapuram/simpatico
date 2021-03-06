<<<<<<< HEAD
#ifndef DDMD_VIRIAL_STRESSTENSOR_AVERAGE_H
#define DDMD_VIRIAL_STRESSTENSOR_AVERAGE_H
=======
#ifndef DDMD_VIRIAL_STRESS_TENSOR_AVERAGE_H
#define DDMD_VIRIAL_STRESS_TENSOR_AVERAGE_H
>>>>>>> 7454788cb955ae6920d2fc65f1ee1fa595a406d8

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/analyzers/Analyzer.h>
#include <ddMd/simulation/Simulation.h>
#include <util/mpi/MpiLoader.h>
#include <util/space/Tensor.h>
#include <util/accumulators/Average.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Periodically write (tensor) StressTensor to file.
   *
   * \ingroup DdMd_Analyzer_Module
   */
   class VirialStressTensorAverage : public Analyzer
   {
   
   public:
   
      /**
      * Constructor.
      *
      * \param simulation parent Simulation object. 
      */
      VirialStressTensorAverage(Simulation& simulation);
   
      /**
      * Destructor.
      */
      virtual ~VirialStressTensorAverage()
      {} 
   
      /**
      * Read dumpPrefix and interval.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
  
      /**
      * Clear nSample counter.
      */
      virtual void clear();
  
      /**
      * Sample virial stress to accumulators
      *
      * \param iStep MD step index
      */
      virtual void sample(long iStep);

      /**
      * Dump configuration to file
<<<<<<< HEAD
      *
      * \param iStep MD step index
=======
>>>>>>> 7454788cb955ae6920d2fc65f1ee1fa595a406d8
      */
      virtual void output();

   private:
 
      /// Output file stream
      std::ofstream outputFile_;
      
      /// Average object to save sxx
      Average sxxAccumulator_;

      /// Average object to save sxx
      Average sxyAccumulator_;

      /// Average object to save sxx
      Average sxzAccumulator_;

      /// Average object to save sxx
      Average syxAccumulator_;

      /// Average object to save sxx
      Average syyAccumulator_;

      /// Average object to save sxx
      Average syzAccumulator_;

      /// Average object to save sxx
      Average szxAccumulator_;

      /// Average object to save sxx
      Average szyAccumulator_;

      /// Average object to save sxx
      Average szzAccumulator_;

      /// Number of samples per block average output
      int nSamplePerBlock_;

      /// Number of configurations dumped thus far (first dump is zero).
      long    nSample_;
   
      /// Has readParam been called?
      long    isInitialized_;
   
   };

}
#endif 
