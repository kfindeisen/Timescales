# Common makefile definitions
# by Krzysztof Findeisen
# Created June 14, 2013
# Last modified June 14, 2013

SHELL := /bin/sh

#---------------------------------------
# Directories
INCDIRS  := 
LIBDIRS  := 

#---------------------------------------
# Compiler settings
CXX       := g++

LANGTYPE  := -std=c++98 -pedantic-errors
WARNINGS  := -Wall -Wextra -Weffc++ -Wdeprecated -Wold-style-cast -Wsign-promo -fdiagnostics-show-option
OPTFLAGS  := -O3 -DNDEBUG
CXXFLAGS  := $(LANGTYPE) $(WARNINGS) $(OPTFLAGS) -Werror -D BOOST_TEST_DYN_LINK
LDFLAGS   := 

#---------------------------------------
# Intermediate builds
AR	:= ar
ARFLAGS	:= rusv

#---------------------------------------
# Other variables
RM      := rm -f
