#ifndef MCMD_MD_SIMULATION_TEST_H
#define MCMD_MD_SIMULATION_TEST_H

#include <test/UnitTestRunner.h>
#include <test/UnitTest.h>

#include <mcMd/mdSimulation/MdSimulation.h>
#include <mcMd/mdSimulation/MdSystem.h>
#include <mcMd/mdIntegrators/MdIntegrator.h>
#include <mcMd/potentials/pair/MdPairPotential.h>
#include <mcMd/potentials/bond/BondPotential.h>
#include <mcMd/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>

using namespace Util;
using namespace McMd;

class MdSimulationTest : public UnitTest
{

public:

   MdSimulationTest()
    : UnitTest(),
      system_(simulation_.system())
   {} 

   virtual void setUp()
   {  simulation_.fileMaster().setRootPrefix(filePrefix()); }

   void testReadParam();
   void testSetZeroVelocities();
   void testSetBoltzmannVelocities();
   void testBuildPairList();
   void testPairEnergy();
   void testAddPairForces();
   void testBondEnergy();
   void testAddBondForces();
   void testCalculateForces();
   void testStep();
   //void testIntegrator();
   void testSimulate();
   void testWriteRestart();
   void testReadRestart();

private:

   MdSimulation simulation_;
   MdSystem&    system_;

};



void MdSimulationTest::testReadParam()
{
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   std::ifstream paramFile;
   openInputFile("in/MdSimulation", paramFile); 
   simulation_.readParam(paramFile);
   simulation_.readCommands();

   try {
      simulation_.isValid();
   } catch (Exception e) {
      TEST_ASSERT(0);
   }

   if (verbose() > 1) {
      std::cout << std::endl;
      simulation_.writeParam(std::cout);
      simulation_.system().writeConfig(std::cout);
   }
}

void MdSimulationTest::testSetZeroVelocities()
{ 
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   std::ifstream paramFile;
   openInputFile("in/MdSimulation", paramFile); 
   simulation_.readParam(paramFile);
   simulation_.readCommands();

   std::cout << std::endl;
   double energy;

   simulation_.simulate(100);
   system_.setZeroVelocities();
   energy = system_.kineticEnergy();
   std::cout << "kinetic energy = " << energy << std::endl;
}

void MdSimulationTest::testSetBoltzmannVelocities()
{ 
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   std::ifstream paramFile;
   openInputFile("in/MdSimulation", paramFile); 
   simulation_.readParam(paramFile);
   simulation_.readCommands();

   std::cout << std::endl;
   //double temperature = 1.0;
   //system_.setBoltzmannVelocities(temperature);
   double energy = system_.kineticEnergy();
   std::cout << "kinetic energy = " << energy << std::endl;
}

void MdSimulationTest::testBuildPairList()
{ 
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   std::ifstream paramFile;
   openInputFile("in/MdSimulation", paramFile); 
   simulation_.readParam(paramFile);
   simulation_.readCommands();

   std::cout << std::endl;

   //double temperature = 1.0;
   //system_.setBoltzmannVelocities(temperature);
   simulation_.simulate(10000);

   try {
      simulation_.isValid();
   } catch (Exception e) {
      std::cout << e.message();
      TEST_ASSERT(0);
   }

}


void MdSimulationTest::testPairEnergy()
{ 
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   std::ifstream paramFile;
   openInputFile("in/MdSimulation", paramFile); 
   simulation_.readParam(paramFile);
   simulation_.readCommands();

   std::cout << std::endl;

   //double temperature = 1.0;
   //system_.setBoltzmannVelocities(temperature);
   simulation_.simulate(1000);

   double energy;
   energy = system_.pairPotential().energy();
   std::cout << energy << std::endl;
}

void MdSimulationTest::testAddPairForces()
{ 
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   std::ifstream paramFile;
   openInputFile("in/MdSimulation", paramFile); 
   simulation_.readParam(paramFile);
   simulation_.readCommands();

   std::cout << std::endl;

   //double temperature = 1.0;
   //system_.setBoltzmannVelocities(temperature);
   simulation_.simulate(1000);
   system_.pairPotential().addForces();
  
}

void MdSimulationTest::testBondEnergy()
{ 
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   std::ifstream paramFile;
   openInputFile("in/MdSimulation", paramFile); 
   simulation_.readParam(paramFile);
   simulation_.readCommands();

   std::cout << std::endl;

   //double temperature = 1.0;
   //system_.setBoltzmannVelocities(temperature);
   simulation_.simulate(1000);

   system_.pairPotential().buildPairList();
   double energy = system_.bondPotential().energy();
   std::cout << energy << std::endl;
}

void MdSimulationTest::testAddBondForces()
{ 
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   std::ifstream paramFile;
   openInputFile("in/MdSimulation", paramFile); 
   simulation_.readParam(paramFile);
   simulation_.readCommands();

   std::cout << std::endl;

   //double temperature = 1.0;
   //system_.setBoltzmannVelocities(temperature);
   simulation_.simulate(1000);

   system_.pairPotential().buildPairList();
   system_.bondPotential().addForces();
}

void MdSimulationTest::testCalculateForces()
{ 
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   std::ifstream paramFile;
   openInputFile("in/MdSimulation", paramFile); 
   simulation_.readParam(paramFile);
   simulation_.readCommands();

   std::cout << std::endl;

   //double temperature = 1.0;
   //simulation_.simulate(1000);
   system_.pairPotential().buildPairList();
   system_.calculateForces();
}

void MdSimulationTest::testStep()
{
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   std::ifstream paramFile;
   openInputFile("in/MdSimulation", paramFile); 
   simulation_.readParam(paramFile);
   simulation_.readCommands();

   std::cout << std::endl;

   double kinetic, potential;
   //double temperature = 1.0;
   simulation_.simulate(100);

   system_.pairPotential().buildPairList();
   system_.calculateForces();
   for (int i=0; i < 10; ++i) {

      kinetic   = system_.kineticEnergy(); 
      potential = system_.potentialEnergy(); 
      std::cout << kinetic << "  " << potential 
                << "  " << kinetic + potential << std::endl;

      system_.mdIntegrator().step();
   }

   kinetic   = system_.kineticEnergy(); 
   potential = system_.potentialEnergy(); 
   std::cout << kinetic << "  " << potential << "  " 
             << kinetic + potential << std::endl;
}

void MdSimulationTest::testSimulate()
{
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   std::ifstream paramFile;
   openInputFile("in/MdSimulation", paramFile); 
   simulation_.readParam(paramFile);
   simulation_.readCommands();

   std::cout << std::endl;

   std::string baseFileName("simulate.0");
   simulation_.writeRestart(baseFileName);

   simulation_.simulate(20);

   baseFileName = "simulate.20";
   simulation_.writeRestart(baseFileName);

}

void MdSimulationTest::testWriteRestart()
{
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   std::ifstream paramFile;
   openInputFile("in/MdSimulation", paramFile); 
   simulation_.readParam(paramFile);
   simulation_.readCommands();

   std::cout << std::endl;

   std::string baseFileName("begin");
   simulation_.writeRestart(baseFileName);

   simulation_.simulate(100000);

   baseFileName = "middle";
   simulation_.writeRestart(baseFileName);

   std::ofstream configFile("middle.cfg");
   simulation_.system().writeConfig(configFile);
   configFile.close();

   bool isContinuation = true;
   simulation_.simulate(100100, isContinuation);

   configFile.open("end.cfg");
   simulation_.system().writeConfig(configFile);
   configFile.close();

}

void MdSimulationTest::testReadRestart()
{
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   std::string baseFileName("middle");
   simulation_.readRestart(baseFileName);

   baseFileName = "middle2";
   simulation_.writeRestart(baseFileName);

   std::ofstream configFile("middle2.cfg");
   simulation_.system().writeConfig(configFile);
   configFile.close();

   bool isContinuation = true;
   simulation_.simulate(100100, isContinuation);

   configFile.open("end2.cfg");
   simulation_.system().writeConfig(configFile);
   configFile.close();
}

#if 0
void MdSimulationTest::testIntegrator()
{
   printMethod(TEST_FUNC);
   std::cout << std::endl;

   double kinetic, potential;
   double temperature = 1.0;
   simulation_.simulate(100);

   system_.pairPotential().buildPairList();
   system_.calculateForces();
   system_.setBoltzmannVelocities(temperature);
   NVEIntegrator integrator(system_);

   for (int i=0; i < 10; ++i) {

      kinetic   = system_.kineticEnergy(); 
      potential = system_.potentialEnergy(); 
      std::cout << kinetic << "  " << potential 
                << "  " << kinetic + potential << std::endl;

      integrator.step();
   }

   kinetic   = system_.kineticEnergy(); 
   potential = system_.potentialEnergy(); 
   std::cout << kinetic << "  " << potential << "  " 
             << kinetic + potential << std::endl;
}
#endif

TEST_BEGIN(MdSimulationTest)
TEST_ADD(MdSimulationTest, testReadParam)
TEST_ADD(MdSimulationTest, testSetZeroVelocities)
TEST_ADD(MdSimulationTest, testSetBoltzmannVelocities)
TEST_ADD(MdSimulationTest, testBuildPairList)
TEST_ADD(MdSimulationTest, testPairEnergy)
TEST_ADD(MdSimulationTest, testAddPairForces)
TEST_ADD(MdSimulationTest, testBondEnergy)
TEST_ADD(MdSimulationTest, testAddBondForces)
TEST_ADD(MdSimulationTest, testCalculateForces)
TEST_ADD(MdSimulationTest, testStep)
   //TEST_ADD(MdSimulationTest, testIntegrator);
TEST_ADD(MdSimulationTest, testSimulate)
TEST_ADD(MdSimulationTest, testWriteRestart)
TEST_ADD(MdSimulationTest, testReadRestart)
TEST_END(MdSimulationTest)

#endif
