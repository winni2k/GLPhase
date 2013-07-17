### Makefile --- for unittesting

## Author: wkretzsch@gmail.com
## Version: $Id: Makefile,v 0.0 2013/06/11 20:46:11 winni Exp $
## Keywords: 
## X-URL: 

# where to find src code to be tested
SRC_DIR= ..

# Points to the root of Google Test, relative to where this file is.
# Remember to tweak this if you move this file.
GTEST_DIR = $(SRC_DIR)/gtest-1.6.0

export CPLUS_INCLUDE_PATH := $(GTEST_DIR)/include:$(SRC_DIR)/snptools.cpp/samtools-0.1.16:$(SRC_DIR)/snptools.cpp/tabix-0.2.5:$(SRC_DIR):$(SRC_DIR)/olivier:$(GTEST_DIR)
export LIBRARY_PATH := $(SRC_DIR)/snptools.cpp/samtools-0.1.16:$(SRC_DIR)/snptools.cpp/tabix-0.2.5
AR= ar

# All tests produced by this Makefile.  Remember to add new tests you
# created to the list.
TESTS = insti_unittest1

test: $(TESTS)

clean:
	rm -f $(TESTS) gtest.a gtest_main.a *.o *.h.gch insti

.PHONY: clean $(TESTS)


# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h


# Builds gtest.a and gtest_main.a.

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
gtest-all.o : $(GTEST_SRCS_)
	$(CC) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CC) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

# Builds a sample test.  A test should link with either gtest.a or
# gtest_main.a, depending on whether it defines its own main()
# function.

insti_unittest1.o : insti_unittest1.cpp $(SRC_DIR)/insti.h $(GTEST_HEADERS)
	$(CC) $(CXXFLAGS) -c insti_unittest1.cpp

insti_unittest1 : $(SRC_DIR)/insti.o $(SRC_DIR)/impute.o $(SRC_DIR)/emcchain.o insti_unittest1.o gtest_main.a $(SRC_DIR)/olivier/utils.o $(SRC_DIR)/relationshipGraph.o
	$(CC) $(CXXFLAGS) -lpthread $^ -o $@ $(LIBS) && ./$@




### Makefile ends here
