#ifndef UTIL_TENSOR_CPP
#define UTIL_TENSOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Tensor.h"
#include "Dimension.h"
#include <util/global.h>

#ifdef UTIL_MPI
#include <util/mpi/MpiStructBuilder.h>
#endif

namespace Util
{

   // Define the Zero Tensor
   const Tensor Tensor::Zero = Tensor(0.0);

   /**
   * Call this function once to guarantee that this file in linked.
   */
   void Tensor::initStatic()
   {
     int nCall = 0;
     if (nCall > 0) UTIL_THROW("Tensor::initStatic called more than once");
     ++nCall;
   }

   // Equality operators
   
   #define UTIL_TENSOR_EPSILON 1.0E-8
   
   bool operator==(const Tensor& t1, const Tensor& t2) 
   {
      for (int i = 0; i < DimensionSq; ++i) {
         if ( fabs(t1.elem_[i] - t2.elem_[i]) > UTIL_TENSOR_EPSILON) {
            return false;
         }
      }
      return true;
   }
   
   bool operator==(const Tensor& t1, const double a2[][Dimension]) 
   {
      for (int i = 0; i < Dimension; ++i) {
         for (int j = 0; j < Dimension; ++j) {
            if ( fabs(t1(i, j) - a2[i][j]) > UTIL_TENSOR_EPSILON) {
               return false;
            }
         }
      }
      return true;
   }
   
   #undef UTIL_TENSOR_EPSILON
   
   bool operator==(const double a1[][Dimension], const Tensor& t2) 
   { return (t2 == a1); }
   
   // Inequality operators
   
   bool operator!=(const Tensor& t1, const Tensor& t2) 
   { return !(t1 == t2); }
   
   bool operator!=(const Tensor& t1, const double a2[][Dimension]) 
   { return !(t1 == a2); }
   
   bool operator!=(const double a1[][Dimension], const Tensor& t2) 
   { return !(t2 == a1); }
   
   /* 
   * Extract a Tensor from an istream.
   */
   std::istream& operator>>(std::istream& in, Tensor &tensor)
   {
      for (int i=0; i < DimensionSq; ++i) {
         in >> tensor.elem_[i];
      }
      return in;
   }
   
   /* 
   * Output a Tensor to an ostream, without line breaks.
   */
   std::ostream& operator<<(std::ostream& out, const Tensor &tensor) 
   {
      for (int i=0; i < DimensionSq; ++i) {
         out.setf(std::ios::scientific);
         out.width(Tensor::Width);
         out.precision(Tensor::Precision);
         out << tensor.elem_[i];
      }
      return out;
   }

   #ifdef UTIL_MPI

   /**
   * Initialize MPI Datatype.
   */
   MPI::Datatype MpiTraits<Tensor>::type = MPI::BYTE;
   bool MpiTraits<Tensor>::hasType = false;

   /**
   * Commit MPI Datatype.
   */
   void Tensor::commitMpiType() 
   {
      if (!MpiTraits<Tensor>::hasType) {
         MpiStructBuilder builder;
         Tensor           tensor;
         int              i, j;
   
         builder.setBase(&tensor);
         for (i = 0; i < Dimension; ++i) {
            for (j = 0; j < Dimension; ++j) {
               builder.addMember(&tensor(i, j), MPI::DOUBLE);
            }
         }
         builder.commit(MpiTraits<Tensor>::type);
         MpiTraits<Tensor>::hasType = true;
      }
   }
   #endif

} 
#endif
