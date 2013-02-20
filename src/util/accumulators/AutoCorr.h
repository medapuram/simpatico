#ifndef UTIL_AUTO_CORR__H
#define UTIL_AUTO_CORR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/accumulators/AutoCorrStage.h>  // base class
#include <util/param/ParamComposite.h>       // base class
#include <util/global.h>
#include <util/format/write.h>

#include <vector>
#include <complex>
using std::complex;

namespace Util
{

   template <typename> class Array;
   
   /**
   * Auto-correlation function for an ensemble of sequences.
   *   
   * This class calculates an autocorrelation function for a ensemble of
   * statistically equivalent sequences x(i) of values of a variable of 
   * type Data. The resulting autocorrelation function is an array of 
   * values of type Product, where C(j) = <x(i-j), x(i)>. Here <A,B> 
   * denotes an inner product of type Product for objects A and B of 
   * type Data.
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
   * \ingroup Accumulators_Module
   */
   template <typename Data, typename Product>
   class AutoCorr : public AutoCorrStage<Data, Product>, public ParamComposite
   {
   
   public:
   
      /// Default constructor.
      AutoCorr();

       /// Default destructor.
      ~AutoCorr();

      /**
      * Reset to empty state.
      */
      void clear();
  
      /**
      * Read parameters, allocate memory and clear history.
      *
      * Reads parameters ensembleCapacity and bufferCapacity, allocates
      * memory, and calls clear().
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Allocate memory, and clear history.
      *
      * Sets parameters ensembleCapacity and bufferCapacity, allocates
      * memory, and calls clear().
      *
      * \param bufferCapacity   maximum number of values in each history
      * \param blockParameter   number of samples per block average
      */
      void setParam(int bufferCapacity, int blockFactor);
   
      void setParam(int bufferCapacity);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);
   
      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
   
      /**
      * Serialize this AutoCorr to/from an archive.
      *
      * \param ar       input or output archive
      * \param version  file version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Sample an array of values and output sampled values (if any).
      *
      * \param values Array of values
      * \param out  output stream to which to write block correlation functions
      */
      void sample(Data value, std::ostream* outFilePtr);

      /**
      * Sample an array of values.
      *
      * If blockFile != 0 and a datafile has been set, this method outputs
      * block correlation functions to the datafile.
      *
      * \param values Array of values
      */
      void sample(Data value);

      /**
      * Output the autocorrelation function
      */
      void output(std::ostream& out);

      /**
      * Return average of all sampled values.
      */
      Data average() const;

      /**
      * Numerical integration of autocorrelation function 
      */
      double corrTime() const;

      /**
      * Return autocorrelation at a given lag time
      * 
      * \param t the lag time
      */
      Product autoCorrelation(int t) const;

      /** 
      * Return maximum number of samples in history for each sequence.
      */
      int bufferCapacity() const;
  
      /**
      * Return blockFactor
      */
      int blockFactor() const;
 
      /** 
      * Return the total number of samples per sequence thus far.
      */
      int nSample() const;

      /**
      * Add pointer to a descendant to an array.
      */
      virtual void registerDescendant(AutoCorrStage<Data, Product>* descendantPtr);

   private:

      /// Array of pointers to descendants.
      std::vector<AutoCorrStage<Data, Product>*>  descendants_;
   
      /// Sum of all previous values of x(t).
      Data sum_;
   
      /// Capacity (# of elements) of corr, nCorr and each buffer.
      int  bufferCapacity_;
      
      int blockFactor_;

      /// Total number of previous values of x(t) per sequence
      int  nSample_;

      /// Private and not implemented to prohibit copying.
      AutoCorr(const AutoCorr& other);
      
      /// Private and not implemented to prohibit assignment.
      AutoCorr& operator = (const AutoCorr& other);
   };

   /*
   * Default constructor.
   */
   template <typename Data, typename Product>
   AutoCorr<Data, Product>::AutoCorr() 
    : AutoCorrStage<Data, Product>(),
      ParamComposite(),
      descendants_(),
      bufferCapacity_(0),
      blockFactor_(0),
      nSample_(0)
   {
      setClassName("AutoCorr"); 
      descendants_.push_back(this);
      setToZero(sum_); 
   }

   /*
   * Destructor.
   */
   template <typename Data, typename Product>
   AutoCorr<Data, Product>::~AutoCorr() 
   {}

   /* 
   * Set accumulator to initial empty state.
   */
   template <typename Data, typename Product>
   void AutoCorr<Data, Product>::clear()
   {   
      setToZero(sum_);
      nSample_ = 0;
      AutoCorrStage<Data, Product>::clear();
   }
 
   /*
   * Read parameters from file.
   */
   template <typename Data, typename Product>
   void AutoCorr<Data, Product>::readParameters(std::istream& in)
   {
      read<int>(in, "bufferCapacity",   bufferCapacity_);
      read<int>(in, "blockFactor",   blockFactor_);
   }
   
   /*
   * Set parameters and initialize.
   */
   template <typename Data, typename Product>
   void AutoCorr<Data, Product>::setParam(int bufferCapacity, int blockFactor)
   {
      bufferCapacity_ = bufferCapacity;
      blockFactor_ = blockFactor;
   }

   /*
   * Set parameters and initialize.
   */
   template <typename Data, typename Product>
   void AutoCorr<Data, Product>::setParam(int bufferCapacity)
   {
      setParam(bufferCapacity, 0);
   }

   /*
   * Load internal state from archive.
   */
   template <typename Data, typename Product>
   void AutoCorr<Data, Product>::loadParameters(Serializable::IArchive &ar)
   {
      loadParameter<int>(ar, "bufferCapacity",   bufferCapacity_);
      loadParameter<int>(ar, "blockFactor",   blockFactor_);
      ar & sum_;
      ar & nSample_;
   }

   /*
   * Save internal state to archive.
   */
   template <typename Data, typename Product>
   void AutoCorr<Data, Product>::save(Serializable::OArchive &ar)
   { ar & *this; }

   /*
   * Serialize this AutoCorr.
   */
   template <typename Data, typename Product>
   template <class Archive>
   void AutoCorr<Data, Product>::serialize(Archive& ar,
                                                const unsigned int version)
   {
      AutoCorrStage<Data, Product>::serialize(ar, version);
      ar & bufferCapacity_;
      ar & blockFactor_;
      ar & sum_;
      ar & nSample_;
   }

   /*
   * Sample a single value from a time sequence.
   */
   template <typename Data, typename Product>
   void AutoCorr<Data, Product>::sample(Data value, std::ostream* outFilePtr)
   {
      AutoCorrStage<Data, Product>::sample(value);

      sum_ += value;

      if (outFilePtr) {
         write<Data>( *outFilePtr, value );
         *outFilePtr << std::endl;
      }
      ++nSample_;
   }

   /*
   * Sample an ensemble of values.
   */
   template <typename Data, typename Product>
   void AutoCorr<Data, Product>::sample(Data value)
   {
      std::ostream* outFilePtr = 0;
      sample(value, outFilePtr);
   }

   /*
   * Final output.
   */
   template <typename Data, typename Product>
   void AutoCorr<Data, Product>::output(std::ostream& outFile) 
   {
      Data    ave;
      Product autocorr;
      Product aveSq;

      ave  = sum_;
      ave /= double(nSample_);
      aveSq = product(ave, ave);

      AutoCorrStage<Data, Product>* ptr = 0;
      int n = descendants_.size();
      for (int i=0; i<n; ++i) {
         ptr = descendants_[i];
         if (i==0) {
            for (int j=0; j < ptr->bufferSize(); ++j) {
               autocorr = ptr->corr(j)/double(ptr->nCorr(j));
               autocorr = autocorr - aveSq;
               outFile << j;
               write<Product>(outFile, autocorr);
               outFile << std::endl;
            }
         } else {
            int minLagIndex = int(bufferCapacity_/blockFactor_);
            for (int j=minLagIndex; j < ptr->bufferSize(); ++j) {
               autocorr = ptr->corr(j)/double(ptr->nCorr(j));
               autocorr = autocorr - aveSq;
               outFile << j * pow((double)blockFactor_, i);
               write<Product>(outFile, autocorr);
               outFile << std::endl;
            }
         }
      }
   }

   /*
   * Return average of sampled values.
   */
   template <typename Data, typename Product>
   Data AutoCorr<Data, Product>::average() const
   {
      Data ave = sum_;
      ave /= double(nSample_);
      return ave;
   }

   /*
   *  Return correlation time in unit of sampling interval  
   */
   template <typename Data, typename Product>
   double AutoCorr<Data, Product>::corrTime() const
   {
      Data    ave;
      Product aveSq, variance, autocorr, sum;

      ave  = sum_;
      ave /= double(nSample_);
      aveSq = product(ave, ave);

      setToZero(sum);

      AutoCorrStage<Data, Product>* ptr = 0;
      int n = descendants_.size();
      long totalTime;
      if (n == 1) {
         totalTime = descendants_[0]->bufferSize() - 1;
      } else {
         totalTime = descendants_[n-1]->bufferSize() - 1;
         totalTime *= pow((double)blockFactor_, n-1);
      } 
      totalTime /= 2;

      long time = 0;
      int i, j;
      for (i=0; i<n; ++i) {
         ptr = descendants_[i];
         if (i==0) {
            variance = ptr->corr(0)/double(ptr->nCorr(0));
            variance = variance - aveSq;
            for (j=1; j < ptr->bufferSize(); ++j) {
               time = j;
               if (time > totalTime)
                  break;
               autocorr = ptr->corr(j)/double(ptr->nCorr(j));
               autocorr = autocorr - aveSq;
               sum += autocorr;
            }
         } else {
            int minLagIndex = int(bufferCapacity_/blockFactor_);
            for (j=minLagIndex; j < ptr->bufferSize(); ++j) {
               time = j*pow((double)blockFactor_, i);
               if (time > totalTime)
                  break;
               autocorr = ptr->corr(j)/double(ptr->nCorr(j));
               autocorr = autocorr - aveSq;
               sum += autocorr;
            }
         }
         if ( j != ptr->bufferSize() ) 
            break;
      }
 
      sum /= variance;
      return sum; 
   }

   /*
   * Return autocorrelation at a given lag time index
   * 
   * \param t the lag time index
   */
   template <typename Data, typename Product>
   Product AutoCorr<Data, Product>::autoCorrelation(int t) const
   {
      Data    ave;
      Product autocorr;
      Product aveSq;

      ave  = sum_;
      ave /= double(nSample_);
      aveSq = product(ave, ave);

      AutoCorrStage<Data, Product>* ptr = 0;
      int n = descendants_.size();
      long totalTime;
      if (n == 1) {
         totalTime = descendants_[n-1]->bufferSize() - 1;
      } else {
         totalTime = descendants_[n-1]->bufferSize() - 1;
         totalTime *= pow((double)blockFactor_, n-1);
      } 
      assert(t < totalTime);

      if (n == 1) {
         ptr = descendants_[0];
         autocorr = ptr->corr(t)/double(ptr->nCorr(t));
         autocorr = autocorr - aveSq;
      } else {
         for (int i=1; i < n; ++i) {
            ptr = descendants_[i];
            int lagMultiplier = pow((double)blockFactor_, i-1);
            int minLagIndex = int(bufferCapacity_/blockFactor_);
            int maxLagIndex = ptr->bufferSize() - 1;
            if (minLagIndex*lagMultiplier <= t <= maxLagIndex*lagMultiplier) {
               int j = int(t/lagMultiplier);
               autocorr = ptr->corr(j)/double(ptr->nCorr(j));
               autocorr = autocorr - aveSq;
               break;
            }
         }
      } 
      return autocorr;
    }
  
   /*
   * Return capacity of history buffer for each sequence.
   */
   template <typename Data, typename Product>
   int AutoCorr<Data, Product>::bufferCapacity() const
   {  return bufferCapacity_; }

   /*
   * Return block factor.
   */
   template <typename Data, typename Product>
   int AutoCorr<Data, Product>::blockFactor() const
   {  return blockFactor_; }

   /*
   * Return number of sampled values.
   */
   template <typename Data, typename Product>
   int AutoCorr<Data, Product>::nSample() const
   {  return nSample_; }

   /*
   * Append pointer to a descendant to descendants_ array.
   */
   template <typename Data, typename Product>
   void AutoCorr<Data, Product>::registerDescendant(AutoCorrStage<Data, Product>* descendantPtr)
   {  descendants_.push_back(descendantPtr); }

}
#endif
