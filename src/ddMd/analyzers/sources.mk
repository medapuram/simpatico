ddMd_analyzers_=\
     ddMd/analyzers/Analyzer.cpp\
     ddMd/analyzers/AnalyzerManager.cpp\
     ddMd/analyzers/AnalyzerFactory.cpp\
     ddMd/analyzers/LogEnergy.cpp\
     ddMd/analyzers/WriteConfig.cpp\
     ddMd/analyzers/OutputEnergy.cpp\
     ddMd/analyzers/OutputPressure.cpp\
     ddMd/analyzers/OutputStressTensor.cpp\
     ddMd/analyzers/VirialStressTensorAverage.cpp\
<<<<<<< HEAD
=======
     ddMd/analyzers/VirialStressTensor.cpp\
>>>>>>> 7454788cb955ae6920d2fc65f1ee1fa595a406d8
     ddMd/analyzers/OutputBoxdim.cpp \
     ddMd/analyzers/OutputTemperature.cpp\
     ddMd/analyzers/OutputPairEnergies.cpp\
     ddMd/analyzers/StructureFactor.cpp\
     ddMd/analyzers/StructureFactorGrid.cpp\
     ddMd/analyzers/VanHove.cpp\
     ddMd/analyzers/OrderParamNucleation.cpp\
     ddMd/analyzers/PairEnergyAverage.cpp\
     ddMd/analyzers/ExternalEnergyAverage.cpp\
     ddMd/analyzers/TrajectoryWriter.cpp\
     ddMd/analyzers/DdMdTrajectoryWriter.cpp\
     ddMd/analyzers/StressAutoCorrelation.cpp

ddMd_analyzers_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_analyzers_))
ddMd_analyzers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(ddMd_analyzers_:.cpp=.o))

