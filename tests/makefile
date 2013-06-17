# Compilation make for timescales test driver
# by Krzysztof Findeisen
# Created June 14, 2013
# Last modified June 17, 2013

include ../makefile.inc

#---------------------------------------
# Select all files
PROJ    := test
SOURCES := driver.cpp unit_stats.cpp unit_lsNormalEdf.cpp unit_FastTable.cpp
OBJS    := $(SOURCES:.cpp=.o)
LIBS    := timescales gsl gslcblas boost_unit_test_framework-mt 

#---------------------------------------
# Primary build option
$(PROJ): $(OBJS)
	@echo "Linking $@ with $(LIBS:%=-l%)"
	@$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(filter %.o,$^) $(DIRS:%=-l%) $(LIBS:%=-l%) $(LIBDIRS:%=-L %) -L .

include ../makefile.common
