
ddMd_sp_processor_= \
    mdCf/processor/Processor.cpp 

# Create lists of source (*.cpp) and object (*.o) files
ddMd_sp_processor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(ddMd_sp_processor_))
ddMd_sp_processor_OBJS=\
     $(addprefix $(OBJ_DIR)/, $(ddMd_sp_processor_:.cpp=.o))

