# compiler variables and flags
CC := g++
# CC := mpic++
CFLAGS := -std=c++11 -g -Wall -O3
# DFLAGS := -DMPI
LIB := -L lib -lyaml-cpp -lboost_system -lboost_filesystem -lfftw3 # -lprofiler
INC := -I include
ifdef FFTW_DIR
    INC += -I $(FFTW_DIR)/include
    LIB += -L $(FFTW_DIR)/lib
endif
TARGET := bin/drift

# build sources and objects lists
SRCDIR := src
BUILDDIR := build
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name '*.$(SRCEXT)')
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

# yaml-cpp install variables
YAML_FILES := lib/libyaml-cpp.a include/yaml-cpp/yaml.h

# boost install variables
BOOST_FILES := lib/libboost_filesystem.a lib/libboost_system.a

# sources for tests
TEST_SRCDIR := test
TEST_BUILDDIR := testbuild
TEST_SOURCES := $(shell find $(TEST_SRCDIR) -name '*.$(SRCEXT)')
TEST_SOURCES += $(filter-out $(SRCDIR)/main.cpp, $(SOURCES))
TEST_OBJECTS := $(addprefix $(TEST_BUILDDIR)/, $(TEST_SOURCES:.$(SRCEXT)=.o))
TEST_LIB_FILES := lib/libgtest.a lib/libgtest_main.a lib/libgmock.a
TEST_LIB_FLAGS := -lgtest -lgtest_main -lgmock -pthread

$(TARGET): $(YAML_FILES) $(BOOST_FILES) $(OBJECTS)
	@echo " Linking..."
	@echo $(SRCDIR)
	@mkdir -p $(dir $@)
	$(CC) $(OBJECTS) -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(DFLAGS) $(INC) -c -o $@ $<

# Install yaml-cpp files in include/ and lib/ under YAML_CPP_INSTALL_PREFIX,
# which is this project's root directory by default
$(YAML_FILES):
	tools/install-yaml-cpp.sh

$(BOOST_FILES):
	tools/install-boost.sh

.PHONY: clean hardclean format

clean:
	@echo " Cleaning...";
	$(RM) -rf $(BUILDDIR)/* $(TEST_BUILDDIR)/* bin/*

hardclean: clean
	$(RM) -rf include/{yaml-cpp,boost,gmock,gtest} lib

format: tools/clang-format-all.sh tools/clang-format
	@echo "Formatting .cpp and .hpp files..."
	$<

tools/clang-format:
	tools/install-clang-format.sh

# Test stuff
.PHONY: test

test: bin/test
	bin/test

bin/test: $(YAML_FILES) $(BOOST_FILES) $(TEST_LIB_FILES) $(TEST_OBJECTS)
	@mkdir -p $(dir $@)
	$(CC) $(TEST_OBJECTS) -o $@ $(LIB) $(TEST_LIB_FLAGS)

$(TEST_LIB_FILES):
	tools/install-googletest.sh

$(TEST_BUILDDIR)/%.o: %.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(DFLAGS) $(INC) -c -o $@ $<
