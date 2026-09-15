include compiler.mk
BUILD_DIR ?= build
MOD_DIR := $(BUILD_DIR)
.DEFAULT_GOAL := check
.PHONY: check clean
check: $(BUILD_DIR)/test_numerics
	./$(BUILD_DIR)/test_numerics
$(BUILD_DIR)/test_numerics: mymats.f90 random.f90 test_numerics.f90 compiler.mk
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) $(MOD_FLAG) -I$(MOD_DIR) -o $@ mymats.f90 random.f90 test_numerics.f90 $(LDLIBS)
clean:
	rm -rf $(BUILD_DIR)
