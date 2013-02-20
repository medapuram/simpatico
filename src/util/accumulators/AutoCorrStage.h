#ifndef UTIL_AUTO_CORR_STAGE_H
#define UTIL_AUTO_CORR_STAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>          // member template
#include <util/containers/RingBuffer.h>      // member template parameter

#include <util/accumulators/setToZero.h>
#include <util/accumulators/product.h>
#include <util/space/Vector.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/global.h>

#include <complex>
using std::complex;

namespace Util
{

   template <typename> class Array;

   /**
   * A hierarchial auto-correlation function for an ensemble of sequences.
   *   
   * This class calculates an autocorrelation function for a ensemble of 
   * statistically equivalent sequences x(i) of values of a variable of 
   * type Data, and also calculates the autocorrelation for sequences of
   * block averages, for multiple block sizes. The algorithm is implemented 
   * by a chain of AutoCorrStage objects, in which each object is 
   * assigned an integer chainId. The "primary" AutoCorrStage object 
   * (chainId=0) in this chain calculates the correlation function for an 
   * ensemble of statistically equivalent "primary" sequences of measured 
   * values that are passed as parameters to its sample method. This object 
   * holds a pointer to an AutoCorrStage with chainId=1 that 
   * calculates the correlation function for an ensemble of statistically 
   * equivalent sequences; A sequence in the ensemble contains values that 
   * are averages of a block of blockFactor consecutive values in the 
   * primary sequence corresponding to the ensemble. The object with 
   * chainId=1 holds a pointer to an object with chainId=2 that calculates 
   * the correlation for an ensemble of sequences of values that are 
   * averages of a block of blockFactor values of sequences processed by 
   * chainId=1, or blockFactor**2 values of the primary sequences. In 
   * general, the object with chainId=n calculates correlations of an 
   * ensemble of sequences of values which are averages of blockFactor**n 
   * values of the primary sequences. The integer parameter blockFactor is 
   * passed to the constructor of the primary AutoCorrStage object as 
   * a parameter, which is set to blockFactor=2 by default. New stages are 
   * added to this chain of objects dynamically as the length of the primary 
   * sequence grows: After the primary AutoCorrStage has sampled m 
   * values, there exists one AutoCorrStage object for each value of 0 <= * 
   * stageId <= n, where n is the largest integer for which blockFactor**n 
   * <= m. 
   *
   * The resulting autocorrelation function at each stage is an array of 
   * values of type Product, where C(j) = <x(i-j), x(i)>. Here <A,B> denotes 
   * an inner product of type Product for objects A and B of type Data.
   *
   * The meaning of <A,B>  for two Data values is defined for various
   * data types by the overloaded functions product(Data, Data) defined
   * in file "product.h" . These functions define a product as an
   * arithmetic product for floating point numbers, and use the 
   * following definitions for complex numbers and Vector objects:
   *
   *     double  product(double, double)   = A*B
   *     complex product(complex, complex) = conjug(A)*B
   *     double  product(Vector,  Vector)  = A.dot(B)
   *
   * The meaning of setting a variable to zero is defined for various
   * types of data by the overloaded functions setToZero(Data&) that 
   * are defined in file setToZero.h.
   *
   * Algorithm:
   *
   * The auto-correlation algorithm is implemented by the sample() method. 
   * Each time the sample method is called with an ensemble of values, the 
   * method updates the buffer of values from which the auto-correlation 
   * can be calculated. It also increments a block sum that is reset to 
   * zero every blockFactor samples. Each AutoCorrStage object can 
   * hold a pointer to a child AutoCorrStage object. When an 
   * AutoCorrStage is instantiated, that pointer is null. After the 
   * first blockFactor samples, however, the sample method creates a new 
   * child object. At this point, and every blockFactor steps thereafter, 
   * the sample method of the parent passes a block average to the sample 
   * method of its child and then resets the block sum to zero. The 
   * resulting chain of objects is extended as needed: After 
   * blockFactor*blockFactor samples have been passed to a parent, so that 
   * blockFactor block average values have been passed to a child object, 
   * the child will create a grandchild. The public interface of the primary 
   * AutoCorrStage object allows all of the objects in the resulting 
   * chain to report error estimates, via the recursive outputError() method, 
   * but does not allow any other form of access to descendants of the 
   * primary AutoCorrStage object.
   *
   * \ingroup Accumulators_Module
   */
   template <typename Data, typename Product>
   class AutoCorrStage
   {
   
   public:
  
      /**
      * Default constructor
      */
      AutoCorrStage();

      /**
      * Destructor.
      *
      * Recursively destroy all children.
      */
      virtual ~AutoCorrStage();

      /**
      * Initialize all accumulators and recursively destroy all children.
      */
      virtual void clear();

      /**
      * Register the creation of a descendant stage.
      *
      * This should be called only by a root stage.
      *
      * \param descendantPtr pointer to a descendant AutoCorrStage.
      */
      virtual void registerDescendant(AutoCorrStage* descendantPtr);
  
      /**
      * Serialize to/from an archive. 
      * 
      * \param ar      archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Sample a value.
      *
      * \param value current value
      */
      void sample(Data value);

      ///\name Accessors
      //@{

      /**
      * Return average of all sampled values.
      */
      Data average() const; 
      
      /**
      * Return calculated correlation 
      * at a lag time (corresponding to lagIndex)
      * 
      * \param lagIndex index of the lag between samples
      */
      Product corr(int lagIndex) const; 

      /**
      * Return number of samples used to calculate correlation 
      * at a lag time (corresponding to lagIndex)
      * 
      * \param lagIndex index of the lag between samples
      */
      int nCorr(int lagIndex) const;

      /** 
      * Return maximum number of samples in history for each sequence.
      */
      int bufferCapacity() const;

      /**
      * Return blockFactor
      */
      int blockFactor() const;

      /**
      * Return size of buffer
      */
      int bufferSize();

      /*
      * Return the number of sampled values.
      */
      long nSample() const;

      /**
      * Return the number of sampled values per block at this stage.
      */
      long stageInterval() const;

      //@}
      
      protected:
      
      /**
      * Does this object have a child AutoCorrStage for block averages?
      */
      bool hasChild() const;
      
      /**
      * Return the child AutoCorrStage by reference.
      */
      AutoCorrStage& child();
      
   private:

      // Ring buffer containing a sequence of stored Data values.
      RingBuffer<Data>  buffer_; 

      // Array in which corr[j] = sum of values of <x(i-j), x(i)>
      DArray<Product>  corr_;
   
      // Array in which nCorr[i] = number of values added to corr[i]
      DArray<int>  nCorr_;
   
      // Sum of all previous values of x(t)
      Data blockSum_;

      // Number of samples used for blockSum evaluation 
      int nBlockSum_;

      // Physical capacity (# of elements) of buffer, corr, and nCorr
      int         bufferCapacity_;
      
      // Physical capacity (# of elements) of buffer, corr, and nCorr
      int         blockFactor_;

      // Total number of previous values of x(t)
      long         nSample_;

      /// Size of buffer
      int  bufferSize_;

      /// Number of measured values per sampled value at this stage.
      long   stageInterval_;
      
      /// Pointer to child stage, if any.
      AutoCorrStage* childPtr_;
      
      /// Pointer to root stage. Null if this is the root stage.
      AutoCorrStage* rootPtr_;
      
      /// Stage index
      int stageId_;
      
      /**
      * Allocate memory and initialize to empty state.
      */
      void allocate();

      /**
      * Is the internal state valid?
      */
      bool isValid();

      /**
      * Constructor for child objects.
      *
      * \param stageInterval number of measured values per sample at this stage
      * \param stageId       integer id for this stage
      * \param rootPtr       pointer to root AverageStage
      * \param blockFactor   ratio of block sizes of subsequent stages
      */
      AutoCorrStage(long stageInterval, int stageId, AutoCorrStage* rootPtr, int bufferCapacity, int blockFactor);

      /**
      * Copy constructor - private and not implemented.
      */
      AutoCorrStage(const AutoCorrStage& other);

      /**
      * Assignment - private and not implemented.
      */
      AutoCorrStage& operator = (const AutoCorrStage& other);

   };
 
   /*
   * Default constructor.
   */
   template <typename Data, typename Product>
   AutoCorrStage<Data, Product>::AutoCorrStage() 
    : buffer_(),
      corr_(),
      nCorr_(),
      blockSum_(),
      nBlockSum_(),
      bufferCapacity_(0),
      blockFactor_(0),
      nSample_(0),
      bufferSize_(0),
      stageInterval_(1),
      childPtr_(0),
      rootPtr_(0),
      stageId_(0)
   {  
      rootPtr_ = this; 
      bufferCapacity_ = rootPtr_->bufferCapacity();
      blockFactor_ = rootPtr_->blockFactor();
      allocate();
   }
   
   /*
   * Default constructor.
   */
   template <typename Data, typename Product>
   AutoCorrStage<Data, Product>::AutoCorrStage(long stageInterval, int stageId, AutoCorrStage* rootPtr, int bufferCapacity, int blockFactor) 
    : buffer_(),
      corr_(),
      nCorr_(),
      blockSum_(),
      nBlockSum_(),
      bufferCapacity_(bufferCapacity),
      blockFactor_(blockFactor),
      nSample_(0),
      bufferSize_(0),
      stageInterval_(stageInterval),
      childPtr_(0),
      rootPtr_(rootPtr),
      stageId_(stageId)
   {  allocate(); }

   template <typename Data, typename Product>
   void AutoCorrStage<Data, Product>::registerDescendant(AutoCorrStage* descendantPtr)
   {}

   /*
   * Destructor.
   */
   template <typename Data, typename Product>
   AutoCorrStage<Data, Product>::~AutoCorrStage()
   {
      if (childPtr_) {
         delete childPtr_;
      }
   }

   /*
   * Set previously allocated to initial empty state.
   */
   template <typename Data, typename Product>
   void AutoCorrStage<Data, Product>::clear()
   {   
      nSample_ = 0;

      if (bufferCapacity_ > 0) {
         nBlockSum_ = 0;
         buffer_.clear();
         setToZero(blockSum_);
         for (int i=0; i < bufferCapacity_; ++i) {
            setToZero(corr_[i]);
            nCorr_[i] = 0;
         }
      }
      if (childPtr_) {
         delete childPtr_;
      }
   }
   
   /*
   * Allocate arrays and CyclicBuffer, and initialize.
   */
   template <typename Data, typename Product>
   void AutoCorrStage<Data, Product>::allocate()
   {  
      buffer_.allocate(bufferCapacity_);
      corr_.allocate(bufferCapacity_);
      nCorr_.allocate(bufferCapacity_);
      clear();
   }
   
   /*
   * Are capacities consistent?
   */
   template <typename Data, typename Product>
   bool AutoCorrStage<Data, Product>::isValid()
   {  
      bool valid = true;
      if (bufferCapacity_ != buffer_.capacity()) valid = false;
      if (bufferCapacity_ != corr_.capacity()) valid = false;
      if (bufferCapacity_ != nCorr_.capacity()) valid = false;
      if (!valid) {
         UTIL_THROW("Invalid AutoCorrStage");
      }
      return valid;
   }
   
   /*
   * Sample a single value from a time sequence.
   */
   template <typename Data, typename Product>
   void AutoCorrStage<Data, Product>::sample(Data value)
   {
      //DArray<Data> blockSumValues;
      //blockSumValues.allocate(ensembleCapacity_);

      buffer_.append(value);
      if (blockFactor_ != 0) {
         blockSum_   += value;
         ++nBlockSum_;
         if (nBlockSum_ == blockFactor_) {
            blockSum_ /= double(blockFactor_);
            if (!childPtr_) {
               long nextStageInterval = stageInterval_*blockFactor_;
               int  nextStageId       = stageId_ + 1;
               childPtr_ = new AutoCorrStage<Data, Product>(nextStageInterval, nextStageId, 
                                                            rootPtr_, bufferCapacity_, blockFactor_);
               rootPtr_->registerDescendant(childPtr_);
            }
            childPtr_->sample(blockSum_);
            setToZero(blockSum_);
            nBlockSum_ = 0;
         }
      }

      if (stageId_ == 0) {
         for (int j = 0; j < buffer_.size(); ++j) {
            corr_[j] += product(buffer_[j], value);
            ++nCorr_[j];
         }
      }
      else {
         int minLagIndex = int(bufferCapacity_/blockFactor_);
         if (buffer_.size() > minLagIndex) {
            for (int j = minLagIndex; j < buffer_.size(); ++j) {
               corr_[j] += product(buffer_[j], value);
               ++nCorr_[j];
            }
         }
      }
      ++nSample_;
   }

   /*
   * Return calculated correlation at a lag time (corresponding to lagIndex)
   */
   template <typename Data, typename Product>
   Product AutoCorrStage<Data, Product>::corr(int lagIndex) const
   {
      if (lagIndex < 0 || lagIndex >= bufferCapacity_)
         UTIL_THROW("lagIndex is out of bounds");

      if (stageId_ == 0) {
         return corr_[lagIndex];
      }
      else {
         int minLagIndex = int(bufferCapacity_/blockFactor_);
         if (lagIndex < minLagIndex)
            UTIL_THROW("lagIndex is less than minimum lag time for this child stage");
         return corr_[lagIndex];
      }
   }
 
   /*
   * Return number of samples used to calculate correlation
   * at a lag time (corresponding to lagIndex)
   */
   template <typename Data, typename Product>
   int AutoCorrStage<Data, Product>::nCorr(int lagIndex) const
   {
      if (lagIndex < 0 || lagIndex >= bufferCapacity_)
          UTIL_THROW("lagIndex is out of bounds");

      if (stageId_ == 0) { 
         return nCorr_[lagIndex];
      }
      else {
         int minLagIndex = int(bufferCapacity_/blockFactor_);
         if (lagIndex < minLagIndex)
            UTIL_THROW("lagIndex is less than minimum lag time for this child stage");
         return nCorr_[lagIndex];
      }
   }

   /*
   * Return capacity of history buffer for each sequence.
   */
   template <typename Data, typename Product>
   int AutoCorrStage<Data, Product>::bufferCapacity() const
   {  return bufferCapacity_; }

   /*
   * Return block factor.
   */
   template <typename Data, typename Product>
   int AutoCorrStage<Data, Product>::blockFactor() const
   {  return blockFactor_; }

   /*
   * Return number of values sampled thus far.
   */
   template <typename Data, typename Product>
   long AutoCorrStage<Data, Product>::nSample() const
   { return nSample_; }

   /*
   * Return size of buffers.
   */
   template <typename Data, typename Product>
   int AutoCorrStage<Data, Product>::bufferSize()
   { 
      bufferSize_ = buffer_.size(); 
      return bufferSize_;
   }

   /*
   * Return the number of measured values per sample at this stage.
   */
   template <typename Data, typename Product>
   long AutoCorrStage<Data, Product>::stageInterval() const
   {  return stageInterval_; }

   /*
   * Serialize this AutoCorrStage.
   */
   template <typename Data, typename Product>
   template <class Archive>
   void AutoCorrStage<Data, Product>::serialize(Archive& ar, 
                                           const unsigned int version)
   {
      ar & bufferCapacity_;
      ar & blockFactor_;
      ar & buffer_;
      ar & corr_;
      ar & nCorr_;
      ar & blockSum_;
      ar & nBlockSum_;
      ar & nSample_;
      ar & stageInterval_;
      isValid();

      int hasChild;
      if (Archive::is_saving()) {
         hasChild = (childPtr_ == 0) ? 0 : 1;
      }
      ar & hasChild;
      if (hasChild) {
         if (Archive::is_loading()) {
            long nextStageInterval = stageInterval_*blockFactor_;
            int  nextStageId       = stageId_ + 1;
            childPtr_ = new AutoCorrStage(nextStageInterval, nextStageId,
                                         rootPtr_, bufferCapacity_, blockFactor_);
            rootPtr_->registerDescendant(childPtr_);
         }
         ar & (*childPtr_);
      } else {
         if (Archive::is_loading()) {
            childPtr_ = 0;
         }
      }


   }

}
#endif
