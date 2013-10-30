# Compilation make for timescales.lib / libtimescales.a
# by Krzysztof Findeisen
# Created March 18, 2010
# Last modified October 29, 2013

include makefile.inc

#-----------------------------------------------------------
# Configurable options

# Do we want to use the fast but high-memory version of lsThreshold, or the slow but low-memory one?
ALGTYPE  := 
#ALGTYPE  := -DSCARGLE_SLOW

#---------------------------------------
# Select all files
PROJ        := timescales
PROJ        := lib$(PROJ).a
SOURCES     := autocorr.cpp dft.cpp freqgen.cpp pairwise.cpp peakfind.cpp scargle.cpp \
	specialfreqs.cpp utils.cpp \
	baddata.cpp badoption.cpp
OBJS        :=     $(SOURCES:.cpp=.o)

#---------------------------------------
# Primary build option
$(PROJ): $(OBJS)
	@echo "Packaging $@..."
	@$(AR) $(ARFLAGS) $@ $^

#---------------------------------------
# Subdirectories
# Can't declare the directories phony directly, or the library will be built every time
.PHONY: cd

examples: cd | $(PROJ)
	@make -C examples --no-print-directory $(MFLAGS)

tests: cd | $(PROJ)
	@make -C tests --no-print-directory $(MFLAGS)

include makefile.common

#---------------------------------------
# Doxygen
.PHONY: doc
doc: doc/
doc/: $(SOURCES) $(DIRS) tests doxygen.cfg
	doxygen doxygen.cfg
	cd doc/latex && make


#---------------------------------------
# Demonstration
.PHONY: example
example: examples | $(PROJ)

#---------------------------------------
# Test cases
.PHONY: unittest
unittest: tests | $(PROJ)

.PHONY: autotest
autotest: $(PROJ) unittest
	@echo "Beginning regression test suite..."
	@echo "Tests started on `date`"
	@cd tests && ./test ; echo "Tests completed on `date`"

#---------------------------------------
# Build program, test suite, and documentation
.PHONY: all
all: $(PROJ) example unittest doc
