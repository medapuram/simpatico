#ifndef UTIL_TEST_CPP
#define UTIL_TEST_CPP

/*
* This program runs all unit tests in the util directory.
*/ 

#ifdef  UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include "accumulators/unit/AccumulatorTestComposite.h"
#include "archives/ArchiveTestComposite.h"
#include "boundary/BoundaryTestComposite.h"
#include "containers/ContainersTestComposite.h"
#include "crystal/CrystalTestComposite.h"
#include "format/FormatTest.h"
#include "param/serial/ParamTestComposite.h"
#include "random/RandomTest.h"
#include "space/SpaceTestComposite.h"

#ifdef  UTIL_MPI
#include "param/mpi/MpiParamTestComposite.h"
#include "mpi/MpiSendRecvTest.h"
#endif

#include <test/CompositeTestRunner.h>

using namespace Util;

TEST_COMPOSITE_BEGIN(UtilNsTestComposite)
addChild(new AccumulatorTestComposite, "accumulators/unit/");
addChild(new ArchiveTestComposite, "archives/");
addChild(new BoundaryTestComposite, "boundary/");
addChild(new ContainersTestComposite, "containers/");
addChild(new CrystalTestComposite, "crystal/");
addChild(new TEST_RUNNER(FormatTest), "format/");
addChild(new ParamTestComposite, "param/serial/");
addChild(new TEST_RUNNER(RandomTest), "random/");
addChild(new SpaceTestComposite, "space/");
#ifdef UTIL_MPI
//addChild(new MpiParamTestComposite, "param/mpi/");
//addChild(new TEST_RUNNER(MpiSendRecvTest), "mpi/");
#endif 
TEST_COMPOSITE_END


int main()
{
   #ifdef UTIL_MPI
   MPI::Init();
   Vector::commitMpiType();
   IntVector::commitMpiType();
   #endif 

   UtilNsTestComposite runner;
   runner.run();


   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif 
}
#endif
