# Compilation make for timescales demo
# by Krzysztof Findeisen
# Created June 14, 2013
# Last modified June 21, 2013

include ../makefile.inc

#---------------------------------------
# Select all files
PROJ    := example
SOURCES := example.cpp
OBJS    := $(SOURCES:.cpp=.o)
LIBS    := timescales gsl

#---------------------------------------
# Primary build option
$(PROJ): $(OBJS)
	@echo "Linking $@ with $(LIBS:%=-l%)"
	@$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(filter %.o,$^) $(DIRS:%=-l%) $(LIBS:%=-l%) $(LIBDIRS:%=-L %) -L .

include ../makefile.common
