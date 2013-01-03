#ifndef DDMD_OUTPUT_BOXDIM_CPP
#define DDMD_OUTPUT_BOXDIM_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "OutputBoxdim.h"
//#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   OutputBoxdim::OutputBoxdim(Simulation& simulation)
    : Diagnostic(simulation),
      nSample_(0),
      isInitialized_(false)
   {}

   /*
   * Read interval and outputFileName.
   */
   void OutputBoxdim::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);

      std::string filename;
      filename  = outputFileName();
      simulation().fileMaster().openOutputFile(filename, outputFile_);
      isInitialized_ = true;
   }

   /*
   * Reset nSample_
   */
   void OutputBoxdim::clear()
   {  nSample_ = 0; }

   /*
   * Dump configuration to file
   */
   void OutputBoxdim::sample(long iStep)
   {
      if (isAtInterval(iStep))  {
         Simulation& sys = simulation();
         sys.computeVirialStress();
         sys.computeKineticStress();
         if (sys.domain().isMaster()) {
            Vector L = sys.boundary().lengths();
            double V = sys.boundary().volume();
            outputFile_ << Int(iStep, 10)
                        << Dbl(L[0], 20)
                        << Dbl(L[1], 20)
                        << Dbl(L[2], 20)
                        << Dbl(V, 20)
                        << std::endl;
         }

         ++nSample_;
      }
   }

}
#endif
