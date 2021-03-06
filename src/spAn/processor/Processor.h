#ifndef SPAN_PROCESSOR_H
#define SPAN_PROCESSOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <spAn/storage/Configuration.h>       // base class
#include <spAn/configIos/ConfigIoFactory.h>   // member 
#include <spAn/analyzers/AnalyzerManager.h>   // member 
#include <util/misc/FileMaster.h>             // member 


namespace SpAn 
{

   class ConfigIo;

   using namespace Util;

   /**
   * A post-processor for analyzing outputs of MD simulations.
   *
   * \ingroup SpAn_Storage_Module
   */
   class Processor : public Configuration
   {

   public:

      /**
      * Constructor
      */
      Processor();

      /**
      * Destructor
      */
      ~Processor();

      /// \name Initialization
      //@{
      
      /**
      * Process command line options.
      *  
      * \param argc number of arguments
      * \param argv array of argument C-strings
      */
      void setOptions(int argc, char * const * argv);

      using ParamComposite::readParam;

      /**
      * Open, read, and close parameter file.
      */
      void readParam(const char* filename);

      /**
      * Read parameters.
      */
      void readParameters(std::istream& in);

      //@}
      /// \name ConfigIo Interface 
      //@{
      
      /**
      * Set ConfigIo style  (creates a ConfigIo).
      */
      void setConfigIo(const std::string& configIoName);

      /**
      * Return the current ConfigIo (create default if necessary).
      */
      ConfigIo& configIo();
   
      /**
      * Read a single configuration file.
      */
      void readConfig(std::ifstream& in);

      /**
      * Open, read and close a configuration file.
      */
      void readConfig(const std::string& filename);
   
      /**
      * Write a single configuration file.
      */
      void writeConfig(std::ofstream& out);

      /**
      * Open, write and close a configuration file.
      */
      void writeConfig(const std::string& filename);

      //@}
      /// \name Trajectory analysis
      //@{

      /**
      * Read and analyze a sequence of numbered configuration files.
      *
      * This function reads and analyzes a sequence of configuration files 
      * that were generated by running a previous simulation. The function 
      * reads files with names of the form inputPrefix() + n for integer 
      * suffixes min <= n <= max. 
      *
      * \param min  integer suffix of first configuration file name
      * \param max  integer suffix of last configuration file name
      * \param fileBaseName root name for dump files (without integer suffix)
      */  
      void analyzeDumps(int min, int max, std::string fileBaseName);

      /**
      * Analyze a trajectory file.
      */
      void analyzeTrajectory(std::string& filename);

      //@}
      /// \name Miscellaneous functions
      //@{
 
      /**
      * Return true if the FileMaster is active.
      */
      bool hasFileMaster() const;

      /**
      * Return FileMaster if active, or throw Exception.
      */
      FileMaster& fileMaster();

      //@}

   private:

      /// Pointer to current ConfigIo object.
      ConfigIo* configIoPtr_;

      /// Factory for generating ConfigIo at run time.
      ConfigIoFactory configIoFactory_;

      /// Manager for analyzers
      AnalyzerManager analyzerManager_;

      /// FileMaster (optionally activated)
      FileMaster fileMaster_;

      /// String identifier for ConfigIo class name
      std::string configIoName_;

   };

}
#endif
