# Compilation make for timescales.lib / libtimescales.a
# by Krzysztof Findeisen
# Created March 18, 2010
# Last modified July 6, 2011

#-----------------------------------------------------------
# Configurable options
CC       = g++
PROJFILE = lib$(PROJ).a
INCDIRS  = c:/krzys/softdev/importlib/include
LIBDIRS  = . /usr/lib c:/krzys/softdev/CPP/lib c:/krzys/softdev/importlib/lib

# Do we want to use the fast but high-memory version of lsThreshold, or the slow but low-memory one?
ALGTYPE  = 
#ALGTYPE  = -DSCARGLE_SLOW

#-----------------------------------------------------------
# Fixed options
PROJ     = timescales
OBJS     = autocorr.o dft.o freqgen.o pairwise.o scargle.o specialfreqs.o utils.o
#TESTOBJS = driver.o unit_autoCorr.o unit_lsNormalEdf.o
TESTOBJS = driver.o unit_stats.o unit_FastTable.o
# Libraries are for test code only. $(PROJ) itself cannot be linked with anything
LIBS	 = gsl gslcblas

# Apparently the Cygwin STL implementation doesn't pass the effective C++ guidelines!
# And Boost.test uses old-style casts
#LANGTYPE = -std=c++98 -Wall -Wextra -Weffc++ -Wdeprecated -Wold-style-cast -Wsign-promo
LANGTYPE = -std=c++98 -Wall -Wextra -Wdeprecated -Wsign-promo
OPTFLAGS = -O3

CFLAGS   = $(LANGTYPE) $(OPTFLAGS) $(ALGTYPE) -Werror

#-----------------------------------------------------------
# Compile instructions
$(PROJFILE): $(OBJS)
	ar rusv $@ $^

%.o: %.cpp
	$(CC) -c $(CFLAGS) -o $@ $^ $(addprefix -I ,$(INCDIRS))
	
# -mno-cygwin linker flag?
example: $(PROJFILE) examples/example.cpp
	$(CC) $(CFLAGS) -o examples/example examples/example.cpp -l$(PROJ) \
		$(addprefix -l,$(LIBS)) $(addprefix -L ,$(LIBDIRS))

test: $(PROJFILE) $(addprefix tests/,$(TESTOBJS))
	$(CC) $(CFLAGS) -o tests/test $(addprefix tests/,$(TESTOBJS)) $(addprefix -L ,$(LIBDIRS)) -l$(PROJ) $(addprefix -l,$(LIBS)) -lboost_unit_test_framework

clean:
	rm -vf *.o *~ examples/*.o examples/*~ tests/*.o tests/*~ core
