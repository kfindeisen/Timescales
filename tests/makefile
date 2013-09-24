# Compilation make for timescales test driver
# by Krzysztof Findeisen
# Created June 14, 2013
# Last modified September 19, 2013

include ../makefile.inc

#---------------------------------------
# Select all files
PROJ    := test
SOURCES := driver.cpp unit_lsNormalEdf.cpp unit_FastTable.cpp
OBJS    := $(SOURCES:.cpp=.o)
LIBS    := kpfutils gsl boost_unit_test_framework-mt 

#---------------------------------------
# Primary build option
$(PROJ): $(OBJS) ../libtimescales.a
	@echo "Linking $@ with ../libtimescales.a $(LIBS:%=-l%)"
	@$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(filter %.o,$^) ../libtimescales.a $(DIRS:%=-l%) $(LIBS:%=-l%) \
		$(LIBDIRS:%=-L %) -L ../../common/ -L .

include ../makefile.common
