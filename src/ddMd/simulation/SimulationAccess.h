#ifndef DDMD_SIMULATION_ACCESS_H
#define DDMD_SIMULATION_ACCESS_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/RArray.h>          // member template
#include <ddMd/chemistry/AtomType.h>         // member template parametr
#include <util/boundary/Boundary.h>          // typedef
#include <ddMd/chemistry/MaskPolicy.h>       // enumeration

namespace Util { 
   class FileMaster;
   class EnergyEnsemble;
   class BoundaryEnsemble;
   class Random;
}

namespace DdMd
{

   class Simulation;
   class AtomStorage;
   template <int N> class GroupStorage;
   class PairPotential;
   class BondPotential;
   #ifdef INTER_ANGLE
   class AnglePotential;
   #endif
   #ifdef INTER_DIHEDRAL
   class DihedralPotential;
   #endif
   #ifdef INTER_EXTERNAL
   class ExternalPotential;
   #endif
   class Domain;
   class Exchanger;
   class AtomType;

   using namespace Util;

   /**
   * Provides access to members of Simulation object.
   *
   * A SimulationAccess holds pointers to the objects owned by a
   * a parent simulation, and values of some variables (such as
   * nAtomType, nBondType, etc.) that are not allowed to change
   * after initialization.
   *
   * It can only be instantiated after all of the objects for
   * which it holds pointers, and after values are set for the
   * variables nAtomType, nBondType, etc.
   *
   * \ingroup DdMd_Simulation_Module
   */
   class SimulationAccess 
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent Simulation object.
      */
      SimulationAccess(Simulation& simulation);

      /**
      * Destructor.
      */
      ~SimulationAccess();

      //@}
      /// \name Accessors (Miscellaneous)
      //@{

     
      /// Get the parent simulation.
      Simulation& simulation();

      /// Get the Boundary.
      Boundary& boundary();
   
      /// Get the AtomStorage.
      AtomStorage& atomStorage();
   
      /// Get the BondStorage.
      GroupStorage<2>& bondStorage();
  
      #ifdef INTER_ANGLE 
      /// Get the AngleStorage.
      GroupStorage<3>& angleStorage();
      #endif
   
      #ifdef INTER_DIHEDRAL
      /// Get the angleStorage.
      GroupStorage<4>& dihedralStorage();
      #endif
   
      /// Get the PairPotential.
      PairPotential& pairPotential();
   
      /// Get the PairPotential.
      BondPotential& bondPotential();
  
      #ifdef INTER_ANGLE 
      /// Get the AnglePotential.
      AnglePotential& anglePotential();
      #endif
   
      #ifdef INTER_DIHEDRAL
      /// Get the DihedralPotential.
      DihedralPotential& dihedralPotential();
      #endif
   
      #ifdef INTER_EXTERNAL
      /// Get the ExternalPotential.
      ExternalPotential& externalPotential();
      #endif
   
      /// Get the EnergyEnsemble.
      EnergyEnsemble& energyEnsemble();

      /// Get the BoundaryEnsemble.
      BoundaryEnsemble& boundaryEnsemble();

      /// Get the Random number generator.
      Random& random();

      /// Get the Domain.
      Domain& domain();

      /// Get the Exchanger.
      Exchanger& exchanger();

      /// Get the FileMaster.
      FileMaster& fileMaster();

      /// Get maximum number of atom types.
      int nAtomType();

      /// Get maximum number of bond types.
      int nBondType();

      #ifdef INTER_ANGLE
      /// Get maximum number of angle types.
      int nAngleType();
      #endif

      #ifdef INTER_DIHEDRAL
      /// Get maximum number of dihedral types.
      int nDihedralType();
      #endif

      #ifdef INTER_EXTERNAL
      /// Does this simulation have an external potential?
      bool hasExternal();
      #endif

      /// Get an AtomType descriptor for atomtype i.
      AtomType& atomType(int i);

      /// Return the value of the mask policy (MaskNone or MaskBonded).
      MaskPolicy maskedPairPolicy() const;

      /// Are forces evaluated by reverse communication (true) or not (false)?
      bool reverseUpdateFlag() const;

      //@}

   private:

      /// Array of AtomType objects for all atoms in a simulation.
      RArray<AtomType> atomTypes_;

      /// Parent Simulation.
      Simulation* simulationPtr_;

      /// Periodic system boundary.
      Boundary* boundaryPtr_;

      /// Container for all atoms and ghosts.
      AtomStorage*  atomStoragePtr_;

      /// Container for bonds.
      GroupStorage<2>*  bondStoragePtr_;

      #ifdef INTER_ANGLE
      /// Container for angles.
      GroupStorage<3>*  angleStoragePtr_;
      #endif

      #ifdef INTER_DIHEDRAL
      /// Container for dihedrals.
      GroupStorage<4>*  dihedralStoragePtr_;
      #endif

      /// Pointer to force/energy evaluator.
      PairPotential* pairPotentialPtr_;

      /// Pointer to covalent bond potential.
      BondPotential* bondPotentialPtr_;

      #ifdef INTER_ANGLE
      /// Pointer to covalent 3-body angle potential.
      AnglePotential* anglePotentialPtr_;
      #endif

      #ifdef INTER_DIHEDRAL
      /// Pointer to covalent 4-body dihedral potential.
      DihedralPotential* dihedralPotentialPtr_;
      #endif

      #ifdef INTER_EXTERNAL
      /// Pointer to external 1-body potential.
      ExternalPotential* externalPotentialPtr_;
      #endif

      /// Pointer to an EnergyEnsemble.
      EnergyEnsemble* energyEnsemblePtr_;

      /// Pointer to an BoundaryEnsemble.
      BoundaryEnsemble* boundaryEnsemblePtr_;

      /// Random number generator.
      Random* randomPtr_;

      /// Processor grid.
      Domain* domainPtr_;

      /// Class to exchange particles, etc.
      Exchanger* exchangerPtr_;

      /// Pointer to a FileMaster.
      FileMaster*  fileMasterPtr_;

      /// Number of distinct atom types.
      int nAtomType_;

      /// Number of distinct bond types.
      int nBondType_;

      #ifdef INTER_ANGLE
      /// Number of distinct angle types.
      int nAngleType_;
      #endif

      #ifdef INTER_DIHEDRAL
      /// Number of distinct dihedral types.
      int nDihedralType_;
      #endif

      #ifdef INTER_EXTERNAL
      /// Does this simulation have an external potential?
      bool hasExternal_;
      #endif

      /**
      * Policy for suppressing pair interactions for some atom pairs.
      *
      * Allowed values of enum MaskPolicy:
      *
      *  - MaskNone:   no masked pairs
      *  - MaskBonded:  mask pair interaction between bonded atoms
      */
      MaskPolicy  maskedPairPolicy_;

      /**
      * Force evaluated using reverse communication (true) or not (false).
      */
      bool reverseUpdateFlag_;

   };

   // Inline method definitions

   inline Simulation& SimulationAccess::simulation()
   {  return *simulationPtr_; }

   inline Boundary& SimulationAccess::boundary()
   {  return *boundaryPtr_; }

   inline AtomStorage& SimulationAccess::atomStorage()
   {  return *atomStoragePtr_; }

   inline GroupStorage<2>& SimulationAccess::bondStorage()
   {  return *bondStoragePtr_; }

   #ifdef INTER_ANGLE
   inline GroupStorage<3>& SimulationAccess::angleStorage()
   {  return *angleStoragePtr_; }
   #endif

   #ifdef INTER_DIHEDRAL
   inline GroupStorage<4>& SimulationAccess::dihedralStorage()
   {  return *dihedralStoragePtr_; }
   #endif

   inline PairPotential& SimulationAccess::pairPotential()
   {  return *pairPotentialPtr_; }

   inline BondPotential& SimulationAccess::bondPotential()
   {  return *bondPotentialPtr_; }

   #ifdef INTER_ANGLE
   inline AnglePotential& SimulationAccess::anglePotential()
   {  return *anglePotentialPtr_; }
   #endif

   #ifdef INTER_DIHEDRAL
   inline DihedralPotential& SimulationAccess::dihedralPotential()
   {  return *dihedralPotentialPtr_; }
   #endif

   #ifdef INTER_EXTERNAL
   inline ExternalPotential& SimulationAccess::externalPotential()
   {  return *externalPotentialPtr_; }
   #endif

   inline EnergyEnsemble& SimulationAccess::energyEnsemble()
   {  return *energyEnsemblePtr_; }

   inline BoundaryEnsemble& SimulationAccess::boundaryEnsemble()
   {  return *boundaryEnsemblePtr_; }

   inline Random& SimulationAccess::random()
   {  return *randomPtr_; }

   inline Domain& SimulationAccess::domain()
   {  return *domainPtr_; }

   inline Exchanger& SimulationAccess::exchanger()
   {  return *exchangerPtr_; }

   inline FileMaster& SimulationAccess::fileMaster()
   {  return *fileMasterPtr_; }

   inline int SimulationAccess::nAtomType()
   {  return nAtomType_; }

   inline int SimulationAccess::nBondType()
   {  return nBondType_; }

   #ifdef INTER_ANGLE
   inline int SimulationAccess::nAngleType()
   {  return nAngleType_; }
   #endif

   #ifdef INTER_DIHEDRAL
   inline int SimulationAccess::nDihedralType()
   {  return nDihedralType_; }
   #endif

   #ifdef INTER_EXTERNAL
   inline bool SimulationAccess::hasExternal()
   {  return hasExternal_; }
   #endif

   inline AtomType& SimulationAccess::atomType(int i)
   {  return atomTypes_[i]; }

   inline MaskPolicy SimulationAccess::maskedPairPolicy() const
   {  return maskedPairPolicy_; }

   inline bool SimulationAccess::reverseUpdateFlag() const
   {  return reverseUpdateFlag_; }

}
#endif
