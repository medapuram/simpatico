#ifndef UTIL_SET_TO_ZERO_H
#define UTIL_SET_TO_ZERO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/space/Vector.h>
#include <util/space/Tensor.h>
#include <complex>

using std::complex;

namespace Util
{

   /**
   * Set a float variable to zero.
   *
   * \param value value to be zeroed.
   */
   inline void setToZero(float& value)
   { value = 0.0; }

   /**
   * Set a double variable to zero.
   *
   * \param value value to be zeroed.
   */
   inline void setToZero(double& value)
   { value = 0.0; }

   /**
   * Set a Vector variable to zero.
   *
   * \param value value to be zeroed.
   */
   inline void setToZero(Vector& value)
   { value.zero(); }

   /**
   * Set a Vector variable to zero.
   *
   * \param value value to be zeroed.
   */
   inline void setToZero(Tensor& value)
   { value.zero(); }

   /**
   * Set a complex<float> variable to zero.
   *
   * \param value value to be zeroed.
   */
   inline void setToZero(complex<float>& value)
   { value = complex<float>(0.0, 0.0); }

   /**
   * Set a complex<double> variable to zero.
   *
   * \param value value to be zeroed.
   */
   inline void setToZero(complex<double>& value)
   { value = complex<double>(0.0, 0.0); }

}
#endif
