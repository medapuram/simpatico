BLD_DIR_REL =../..
include $(BLD_DIR_REL)/config.mk
include $(BLD_DIR)/spAn/config.mk
include $(BLD_DIR)/inter/config.mk
include $(BLD_DIR)/util/config.mk
include $(SRC_DIR)/spAn/patterns.mk
include $(SRC_DIR)/inter/sources.mk
include $(SRC_DIR)/util/sources.mk
include $(SRC_DIR)/spAn/sources.mk
include $(SRC_DIR)/spAn/tests/sources.mk


TEST=spAn/tests/Test

all: $(spAn_tests_OBJS)

run: $(spAn_tests_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/spAn/tests/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in spAn/tests/log" > count
	@cat count

clean:
	rm -f $(spAn_tests_OBJS) $(spAn_tests_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST)
	cd chemistry; $(MAKE) clean
	cd configIos; $(MAKE) clean
	cd neighbor; $(MAKE) clean
	cd storage; $(MAKE) clean
	cd processor; $(MAKE) clean
	$(MAKE) clean-outputs

clean-outputs:
	rm -f log count
	#cd chemistry; $(MAKE) clean-outputs
        
-include $(spAn_tests_OBJS:.o=.d)
-include $(spAn_OBJS:.o=.d)
-include $(inter_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)


