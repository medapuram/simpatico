#ifndef UTIL_SERIALIZABLE_H
#define UTIL_SERIALIZABLE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace Util
{

   class BinaryFileOArchive;
   class BinaryFileIArchive;

   /**
   * Abstract class for serializable objects.
   *
   * This class defines an interface for serialization of objects.
   * The save() method saves the internal state of an object to an 
   * archive, and the load() method loads the state from an archive.
   *
   * The type of archive to be used is specified by the OArchiveType 
   * and IArchiveType typedefs.  The two concrete classes that are 
   * referred to by these typedefs should be forward declared in 
   * this file, and the header files for these two classes must be
   * included in the file Serializable_includes.h. The file
   * Serializable_includes.h should be included in source files
   * that implement that load and save methods for subclasses.
   */
   class Serializable
   {

   public:

      /**
      * Type of output archive used by save method.
      */
      typedef BinaryFileOArchive OArchiveType;

      /**
      * Type of input archive used by load method.
      */
      typedef BinaryFileIArchive IArchiveType;

      /**
      * Destructor.
      */
      virtual ~Serializable(){};

      /**
      * Save to an archive.
      *
      * \param ar binary saving (output) archive.
      */
      virtual void save(OArchiveType& ar) = 0;

      /**
      * Load from an archive.

      * \param ar binary loading (input) archive.
      */
      virtual void load(IArchiveType& ar) = 0;

   };

}
#endif
