
##
## This file builds assignment 2
##

##
## Some macros (makefile 'variables') are defined by default
## by make(1).  The ones relevant to C compiling are:
##  CC     : the name of the compiler (defaults to system compiler)
##  CFLAGS : the flags (options) to pass to the compiler
##
## There are similar macros for other languages, for example C++ uses
##  CXX      : the name of the c++ compiler (defaults to system compiler)
##  CXXFLAGS : the flags (options) to pass to the c++ compiler


## explicitly add debugger support to each file compiled,
## and turn on all warnings.  If your compiler is surprised by your
## code, you should be too.
CFLAGS = -g -Wall

## uncomment/change this next line if you need to use a non-default compiler
#CC = cc


##
## We can define variables for values we will use repeatedly below
##

## define the executables we want to build
LOEXE = llloadonly
HOEXE = llheadonly
HTEXE = llheadtail
ADEXE = arraydouble

## Define the set of object files we need to build each executable.
## If you write more files, be sure to add them in here
LOOBJS		= llloadonly_main.o fasta_read.o
HOOBJS		= llheadonly_main.o LLvNode.o fasta_read.o
HTOBJS		= llheadtail_main.o LLvNode.o fasta_read.o
ADOBJS		= arraydouble_main.o fasta_read.o


##
## TARGETS: below here we describe the target dependencies and rules
##
all: $(LOEXE) $(HOEXE) $(HTEXE) $(ADEXE)

$(HOEXE): $(HOOBJS)
	$(CC) $(CFLAGS) -o $(HOEXE) $(HOOBJS)

$(LOEXE): $(LOOBJS)
	$(CC) $(CFLAGS) -o $(LOEXE) $(LOOBJS)

$(HTEXE): $(HTOBJS)
	$(CC) $(CFLAGS) -o $(HTEXE) $(HTOBJS)

$(ADEXE): $(ADOBJS)
	$(CC) $(CFLAGS) -o $(ADEXE) $(ADOBJS)

## convenience target to remove the results of a build
clean :
	- rm -f $(LOOBJS) $(LOEXE)
	- rm -f $(HOOBJS) $(HOEXE)
	- rm -f $(HTOBJS) $(HTEXE)
	- rm -f $(ADOBJS) $(ADEXE)

