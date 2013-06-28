#ifndef DDMD_CONFIG_IO_H
#define DDMD_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>           // base class
#include <ddMd/communicate/AtomDistributor.h>    // member 
#include <ddMd/communicate/AtomCollector.h>      // member 
#include <ddMd/communicate/GroupDistributor.h>   // member 
#include <ddMd/communicate/GroupCollector.h>     // member 
#include <util/boundary/Boundary.h>              // typedef

#include <util/containers/DArray.h>              // member

#include <ddMd/chemistry/MaskPolicy.h>

namespace DdMd
{

   class Simulation;
   class Domain;
   class AtomStorage;
   class BondStorage;
   class AngleStorage;
   class DihedralStorage;
   class Buffer;

   using namespace Util;

   /**
   * Abstract reader/writer for configuration files.
   *
   * Each concrete subclass of ConfigIo implements a specific file format
   * by implementing the pure virtual readConfig and writeConfig methods. 
   *
   * \ingroup DdMd_ConfigIo_Module
   */
   class ConfigIo  : public ParamComposite
   {

   public:

      /**
      * Default constructor.
      */
      ConfigIo();

      /**
      * Constructor.
      *
      * \param simulation parent Simulation object.
      */
      ConfigIo(Simulation& simulation);

      /**
      * Associate with related objects.
      *
      * Required iff instantiated with default constructor.
      */
      void associate(Domain& domain, Boundary& boundary,
                     AtomStorage& atomStorage,
                     BondStorage& bondStorage,
                     #ifdef INTER_ANGLE
                     AngleStorage& angleStorage,
                     #endif
                     #ifdef INTER_DIHEDRAL
                     DihedralStorage& dihedralStorage,
                     #endif
                     Buffer& buffer);

      /**
      * Set cache sizes and allocate memory.
      *
      * \param atomCacheCapacity size of internal atom cache. 
      * \param bondCacheCapacity size of internal bond cache. 
      * \param angleCacheCapacity size of internal angle cache. 
      * \param dihedralCacheCapacity size of internal dihedral cache. 
      */
      virtual void initialize(int atomCacheCapacity = 100,
                              int bondCacheCapacity = 100
                              #ifdef INTER_ANGLE
                              , int angleCacheCapacity = 100
                              #endif
                              #ifdef INTER_DIHEDRAL                             
                              , int dihedralCacheCapacity = 100
                              #endif
                              );

      /**
      * Read cache size and allocate memory.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load internal state of configIo (not configuration) from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void load(Serializable::IArchive &ar);

      /**
      * Save internal state of ConfigIo (not configuration) to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Read a configuration file.
      *
      * This method reads a configuration file from the master processor, 
      * and distributes atoms and groups among the processors. It should
      * be implemented using an AtomDistributor object. When Atoms are 
      * added to the AtomDistributor on the master processor, the atomic
      * coordinates should be in generalized / scaled coordinates, such 
      * that 0 <= position[i] < 1.0 for i =0,..,2. Atomic positions that
      * lie slightly outside the domain [0.0 ,1.0]  will automatically 
      * be shifted into this domain by the AtomDistributor::addAtom() 
      * method before being assigned to processors or distributed.
      *
      * Upon entry (preconditions):
      *
      *  - The configuration file must be open for reading.
      *
      *  - The AtomStorage must be set for generalized / scaled coordinates.
      *
      * Upon return (postconditions):
      *
      *   - Each processor owns all atoms in its domain, and no others.
      *     Each atom is owned by one and only one processor.
      *
      *   - Each atom position is in the primary simulation cell.
      *
      *   - Atom coordinates are in scaled form, between 0.0 and 1.0.
      *
      *   - Each processor owns every group that contains one or more of 
      *     the atoms that it owns, and no empty groups. A group may be
      *     ``owned" by more than one processor.
      *
      *   - All Atom Mask data is set correctly for specified maskPolicy.
      *
      *   - There are no ghost atoms on any processor.
      *
      * \param file input file stream (must be open on master).
      * \param maskPolicy MaskPolicy to be used in setting atom masks
      */
      virtual void readConfig(std::ifstream& file, MaskPolicy maskPolicy) = 0;

      /**
      * Write configuration file.
      *
      * This function writes a file on the master processor by collecting 
      * atom and group data from all processors. Implementations must 
      * work correctly when atom positions stored in either Cartesian 
      * or scaled / generalized coordinates on entry and exit, by testing 
      * AtomStorage::isCartesian(), and converting coordinate systems as 
      * needed before writing each position file. Atomic positions may be
      * written in Cartesian or generalized coordinates, depending on the
      * file format, but the same convention should be used in writeConfig() 
      * and readConfig() in each class. Some file formats may allow atom and 
      * groups to be written in arbitrary order. Atomic positions may be 
      * written "as is", and may lie slightly outside the primary cell, 
      * because they will be shifted back by the AtomDistributor used to 
      * implement readConfig(). This function may not modify actual Atom
      * positions or change the AtomStorage isCartesian() flag.
      *
      * Usage:
      *
      *   - When used within Simulation::readCommands() to implement the 
      *     WRITE_CONFIG command, atomic coordinates will be scaled / 
      *     generalized on both entry and exit. 
      *
      *   - When used within a Diagnostic, to dump configuration files, 
      *     atomic coordinates will be Cartesian on entry and exit.
      *
      * \param file output file stream (must be open on master processor).
      */
      virtual void writeConfig(std::ofstream& file) = 0;

   protected:

      /**
      * Set Mask data on all atoms.
      */
      void setAtomMasks();

      /**
      * Get the Domain by reference.
      */
      Domain& domain();

      /**
      * Get Boundary by reference.
      */
      Boundary& boundary();
   
      /**
      * Get AtomStorage by reference.
      */
      AtomStorage& atomStorage();
   
      /**
      * Get the AtomDistributor by reference.
      */
      AtomDistributor& atomDistributor();

      /**
      * Get the AtomCollector by reference.
      */
      AtomCollector& atomCollector();

      /**
      * Get BondStorage by reference.
      */
      BondStorage& bondStorage();
  
      /**
      * Get the AtomDistributor by reference.
      */
      GroupDistributor<2>& bondDistributor();

      /**
      * Get the bond collector by reference.
      */
      GroupCollector<2>& bondCollector();

      #ifdef INTER_ANGLE
      /**
      * Get AngleStorage by reference.
      */
      AngleStorage& angleStorage();

      /**
      * Get the angle distributor by reference.
      */
      GroupDistributor<3>& angleDistributor();

      /**
      * Get the angle collector by reference.
      */
      GroupCollector<3>& angleCollector();
      #endif

      #ifdef INTER_DIHEDRAL
      /**
      * Get DihedralStorage by reference.
      */
      DihedralStorage& dihedralStorage();

      /**
      * Get the dihedral distributor by reference.
      */
      GroupDistributor<4>& dihedralDistributor();

      /**
      * Get the dihedral collector by reference.
      */
      GroupCollector<4>& dihedralCollector();
      #endif

   private:

      // Distributors and collectors
      AtomDistributor atomDistributor_;
      AtomCollector atomCollector_;
      GroupDistributor<2> bondDistributor_;
      GroupCollector<2> bondCollector_;
      #ifdef INTER_ANGLE
      GroupDistributor<3> angleDistributor_;
      GroupCollector<3> angleCollector_;
      #endif
      #ifdef INTER_DIHEDRAL
      GroupDistributor<4> dihedralDistributor_;
      GroupCollector<4> dihedralCollector_;
      #endif

      // Pointers to associated objects.
      Domain* domainPtr_;
      Boundary* boundaryPtr_;
      AtomStorage* atomStoragePtr_;
      BondStorage* bondStoragePtr_;
      #ifdef INTER_ANGLE
      AngleStorage* angleStoragePtr_;
      #endif
      #ifdef INTER_DIHEDRAL
      DihedralStorage* dihedralStoragePtr_;
      #endif

      // Cache capacities
      int  atomCacheCapacity_;
      int  bondCacheCapacity_;
      #ifdef INTER_ANGLE
      int  angleCacheCapacity_;
      #endif
      #ifdef INTER_DIHEDRAL
      int  dihedralCacheCapacity_;
      #endif

      /**
      * Read Group<N> objects from file. 
      */
      template <int N>
      int readGroups(std::ifstream& file, 
                     const char* sectionLabel, const char* nGroupLabel,
                     GroupDistributor<N>& distributor);

      /**
      * Write Group<N> objects to file. 
      */
      template <int N>
      int writeGroups(std::ofstream& file, 
                      const char* sectionLabel, const char* nGroupLabel,
                      GroupStorage<N>& storage, GroupCollector<N>& collector);

   };

   // Inline method definitions

   inline Domain& ConfigIo::domain()
   {  return *domainPtr_; }

   inline Boundary& ConfigIo::boundary()
   {  return *boundaryPtr_; }

   inline AtomStorage& ConfigIo::atomStorage()
   {  return *atomStoragePtr_; }

   inline AtomDistributor& ConfigIo::atomDistributor()
   {  return atomDistributor_; }

   inline AtomCollector& ConfigIo::atomCollector()
   {  return atomCollector_; }

   inline BondStorage& ConfigIo::bondStorage()
   {  return *bondStoragePtr_; }

   inline GroupDistributor<2>& ConfigIo::bondDistributor()
   {  return bondDistributor_; }

   inline GroupCollector<2>& ConfigIo::bondCollector()
   {  return bondCollector_; }

   #ifdef INTER_ANGLE
   inline AngleStorage& ConfigIo::angleStorage()
   {  return *angleStoragePtr_; }

   inline GroupDistributor<3>& ConfigIo::angleDistributor()
   {  return angleDistributor_; }

   inline GroupCollector<3>& ConfigIo::angleCollector()
   {  return angleCollector_; }
   #endif

   #ifdef INTER_DIHEDRAL
   inline DihedralStorage& ConfigIo::dihedralStorage()
   {  return *dihedralStoragePtr_; }

   inline GroupDistributor<4>& ConfigIo::dihedralDistributor()
   {  return dihedralDistributor_; }

   inline GroupCollector<4>& ConfigIo::dihedralCollector()
   {  return dihedralCollector_; }
   #endif

}
#endif
