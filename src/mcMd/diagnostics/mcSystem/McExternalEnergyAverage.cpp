#ifdef INTER_EXTERNAL
#ifndef MCMD_MC_EXTERNAL_ENERGY_AVERAGE_CPP
#define MCMD_MC_EXTERNAL_ENERGY_AVERAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McExternalEnergyAverage.h"        // class header

#include <util/misc/FileMaster.h>  
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/potentials/external/ExternalPotential.h>


#include <cstdio> 

namespace McMd
{

   using namespace Util;

   /* 
   * Constructor.
   */
   McExternalEnergyAverage::McExternalEnergyAverage(McSystem& system)
    : AverageDiagnostic<McSystem>(system)
   {  setClassName("McExternalEnergyAverage"); }
 
   /*
   * Evaluate external energy, and add to ensemble. 
   */
   void McExternalEnergyAverage::sample(long iStep) 
   {
      if (!isAtInterval(iStep)) return;
      int n = 0;
      double energy = 0.0;
      System::MoleculeIterator molIter;
      Molecule::AtomIterator atomIter;
      for (int iSpec=0; iSpec < system().simulation().nSpecies(); ++iSpec){
         for (system().begin(iSpec, molIter); molIter.notEnd(); ++molIter){
             for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
                 energy += system().externalPotential().energy(atomIter->position(), atomIter->typeId());
                 if (n == 0){
                    std::cout << "energy for 1st atom is " <<  system().externalPotential().energy(atomIter->position(), atomIter->typeId()) << std::endl;
                    std::cout << "position for 1st atom is " << atomIter->position() << std::endl;
                 }
             ++n;
             }
         }
      }
      accumulator_.sample(energy, outputFile_);
   }

}
#endif
#endif 
