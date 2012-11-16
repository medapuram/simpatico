#ifndef INTER_ORDERING_EXTERNAL_H
#define INTER_ORDERING_EXTERNAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, Jian Qin and David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/boundary/Boundary.h>
#include <util/space/Dimension.h>
#include <util/param/ParamComposite.h>
#include <util/global.h>
#include <cmath>

namespace Inter
{

   using namespace Util;

   /**
   * A clipped cosine potential that induces lamellar ordering
   * along the direction specified by perpDirection_.
   *    perpDirection_ = 0: x direction
   *                   = 1: y direction
   *                   = 2: z direction
   *
   *                                                  /                   /                   z   \ \
   *  u = prefactor[atomType] externalParameter tanh | clipParameter cos | 2  pi periodicity ---   | |
   *                                                  \                   \                   Lz  / / 
   *
   * Prefactor (which depends on the atomType), externalParameter, interfaceWidth (relative to the box length 
   * along the direction perpendicular to lamellae) and periodicity are given as inputs in the parameter file. 
   * ClipParameter is the inverse of 2*pi*periodicity*interfaceWidth. 
   *
   * \ingroup Inter_External_Module
   */
   class OrderingExternal : public ParamComposite 
   {
   
   public:
   
      /**
      * Default constructor.
      */
      OrderingExternal();

      /**
      * Copy constructor.
      */
      OrderingExternal(const OrderingExternal& other);

      /**
      * Assignment.
      */
      OrderingExternal& operator = (const OrderingExternal& other);

      /**  
      * Set nAtomType value.
      *
      * \param nAtomType number of atom types.
      */
      void setNAtomType(int nAtomType);

      /**
      * Sets external parameter
      *
      * \param externalParameter external parameter of system
      */
      void setExternalParameter(double externalParameter);

      /**
      * Set pointer to Boundary.
      *
      * \param boundary Boundary object (used to calculate length along perpendicular direction).
      */
      void setBoundary(Boundary &boundary);

      /**
      * Read potential parameters, and initialize other variables.
      *
      * \pre nAtomType must have been set, by calling setNAtomType().
      * \pre Boundary must have been set, by calling setBoundary().
      *
      * \param in input stream 
      */
      void readParameters(std::istream &in);

      /**
      * Returns external parameter
      *
      * \return external parameter
      */
      double externalParameter() const;

      /**
      * Returns external potential energy of a particle of type i.
      *
      * \param d  component of position along the perpendicular direction
      * \param i  type of particle (prefactor depends on atomtype)
      * \return   external potential energy
      */
      //double energy(double d, int i) const;
 
      /**
      * Returns external potential energy of a single particle. 
      *
      * \param position atomic position Vector
      * \param i        atom type.
      * \return     external potential energy
      */
      double energy(const Vector& position, int i) const;

      /**
      * Returns magnitude of the external force.
      *
      * \param d    component of position along the perpendicular direction
      * \param type atom type id (not used)
      * \return    force scalar
      */
      //double forceScalar(double d, int type) const;
 
      /**
      * Returns force caused by the external potential.
      *
      * \param position  atom position
      * \param type      atom type id
      * \param force     force on the atom (on output)
      */
      void getForce(const Vector& position, int type, Vector& force) const;
 
      /**
      * Return name string "OrderingExternal".
      */
      std::string className() const;
 
   private:
   
      /// Maximum allowed value for nAtomType (# of particle types).
      static const int MaxAtomType = 2;

      /// Prefactor array ofsize nAtomType.
      DArray<double> prefactor_;

      /// External parameter.
      double externalParameter_;

      /// Number of reciprocal lattice vectors
      int  nWaveVectors_;

      /// Array of Miller index IntVectors for the reciprocal lattice vectors.
      DArray<IntVector>  waveIntVectors_;

      /// Array of floating point reciprocal lattice vectors.
      DArray<Vector>  waveVectors_;

      /// Pointer to associated Boundary object.
      Boundary *boundaryPtr_;
   
      /// Number of possible atom types.
      int    nAtomType_; 

      /// Are all parameters and pointers initialized?
      bool  isInitialized_;

   };
  
   // inline methods 

   /* 
   * Calculate external potential energy for a single atom.
   */
   inline double OrderingExternal::energy(const Vector& position, int type) const
   {
      Vector cellLengths = boundaryPtr_->lengths();
      Vector reducedPosition = position;
      for (int i = 0; i < Dimension; ++i) {
         reducedPosition[i] /= cellLengths[i];
      }

      DArray<double> arg;
      arg.allocate(nWaveVectors_);

      double cosine = 0.0;

      for (int i = 0; i < nWaveVectors_; ++i) {
         arg[i] = reducedPosition.dot(waveVectors_[i]);
         cosine += cos(arg[i]);
      }

      return prefactor_[type]*externalParameter_*tanh(cosine);
   }

   /* 
   * Calculate external force for a single atom.
   */
   inline 
   void OrderingExternal::getForce(const Vector& position, int type, 
                                     Vector& force) const
   {
      DArray<double> arg;
      arg.allocate(nWaveVectors_);

      DArray<double> sine;
      sine.allocate(nWaveVectors_);

      Vector deriv;
      deriv.zero();

      DMatrix<double> derivWave;
      derivWave.allocate(nWaveVectors_, Dimension);

      double cosine = 0.0;
      int i,j;

      Vector dimensions = boundaryPtr_->lengths();
      Vector reducedPosition = position;
      for (j = 0; j < Dimension; ++j) {
         reducedPosition[j] /= dimensions[j];
      }

      for (i = 0; i < nWaveVectors_; ++i) {
         arg[i] = reducedPosition.dot(waveVectors_[i]);
         sine[i] = sin(arg[i]);

         for (j = 0; j < Dimension; ++j) {
            derivWave(i,j) = -sine[i]*waveVectors_[i][j];
            deriv[j] += derivWave(i,j);
         }

         cosine += cos(arg[i]);
      }

      double tanH = tanh(cosine);
      double sechSq = (1.0 - tanH*tanH);
      force.zero();

      for (j = 0; j < Dimension; ++j) {
         force[j] += prefactor_[type]*externalParameter_*sechSq;
         force[j] *= deriv[j];
      }
   }
 
}
#endif
