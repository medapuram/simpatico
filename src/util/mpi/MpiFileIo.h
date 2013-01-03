#ifndef UTIL_MPI_FILE_IO_H
#define UTIL_MPI_FILE_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include <iostream>
#include <string>

namespace Util
{

   /**
   * Identifies whether this processor may do file I/O.
   *
   * The return value of isIoProcessor() indicates whether this processor
   * can read and write to file. If the the class is compiled with UTIL_MPI
   * not defined, then isIoProcessor() always returns true. If the class is
   * compiled with UTIL_MPI defined, then this method returns true if either:
   * (1) A communicator has been set  * and this processor has rank 0 within 
   * that communicator, or (2) No communicator has been set.
   *
   * When compiled with UTIL_MPI defined, an MpiFileIo object has a pointer 
   * to an MPI communiciator, and provides methods to set and unset (nullify) the
   * associated communicator.
   */
   class MpiFileIo
   {

   public:

      /**
      * Constructor.
      */
      MpiFileIo();

      /**
      * Copy constructor.
      */
      MpiFileIo(const MpiFileIo& other);

      /**
      * Can this processor do file I/O ?
      */
      bool isIoProcessor() const;

      #ifdef UTIL_MPI
      /**
      * Set the  communicator.
      */
      void setCommunicator(MPI::Intracomm& communicator);

      /**
      * Clear (nullify) the  communicator.
      */
      void clearCommunicator();

      /**
      * Does this object have a  communicator?
      */
      bool hasCommunicator() const;

      /**
      * Get the  communicator.
      */
      MPI::Intracomm& communicator() const;
      #endif

   private:

      /// Should this object read and write?
      bool isIoProcessor_;

      #ifdef UTIL_MPI
      /// Pointer to the  communicator.
      MPI::Intracomm* communicatorPtr_;
      #endif

   };

   // Inline methods

   /*
   * Should this processor do file I/O ?
   */
   inline bool MpiFileIo::isIoProcessor() const
   {  return isIoProcessor_; }

   #ifdef UTIL_MPI

   inline bool MpiFileIo::hasCommunicator() const
   {  return (communicatorPtr_ != 0); }

   /*
   * Get the  communicator.
   */
   inline MPI::Intracomm& MpiFileIo::communicator() const
   {
      assert(communicatorPtr_);  
      return *communicatorPtr_; 
   }
   #endif
}

#endif
