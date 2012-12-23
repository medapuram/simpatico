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
#include <util/space/Vector.h>
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

      /// Interface widths array ofsize nWaveVectors
      DArray<double> interfaceWidths_;

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
      const Vector cellLengths = boundaryPtr_->lengths();

      double cosine = 0.0;

      for (int i = 0; i < nWaveVectors_; ++i) {
         Vector q;
         q[0] = 2.0*M_PI*waveIntVectors_[i][0]/cellLengths[0];
         q[1] = 2.0*M_PI*waveIntVectors_[i][1]/cellLengths[1]; 
         q[2] = 2.0*M_PI*waveIntVectors_[i][2]/cellLengths[2];
         double arg, clipParameter;
         arg = position[0]*q[0] + position[1]*q[1] + position[2]*q[2];
         clipParameter = 1.0/( interfaceWidths_[i]*( q[0]*cellLengths[0] + q[1]*cellLengths[1] + q[2]*cellLengths[2] ) );
         cosine += clipParameter*cos(arg);
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
      const Vector cellLengths = boundaryPtr_->lengths();

      double cosine = 0.0;
      Vector deriv;
      deriv.zero();
      for (int i = 0; i < nWaveVectors_; ++i) {
         Vector q;
         q[0] = 2.0*M_PI*waveIntVectors_[i][0]/cellLengths[0];
         q[1] = 2.0*M_PI*waveIntVectors_[i][1]/cellLengths[1]; 
         q[2] = 2.0*M_PI*waveIntVectors_[i][2]/cellLengths[2];
         double arg, sine, clipParameter;
         arg = position[0]*q[0] + position[1]*q[1] + position[2]*q[2];
         clipParameter = 1.0/( interfaceWidths_[i]*( q[0]*cellLengths[0] + q[1]*cellLengths[1] + q[2]*cellLengths[2] ) );
         cosine += clipParameter*cos(arg);
         sine = clipParameter*sin(arg);
         q *= sine;
         deriv -= q;
         cosine += clipParameter*cos(arg);
      }
      double tanH = tanh(cosine);
      double sechSq = (1.0 - tanH*tanH);
      double f = prefactor_[type]*externalParameter_*sechSq;
      deriv *= f;
      force = deriv;
   }
 
}
#endif
