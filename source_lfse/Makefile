##
###################################################################
##                                                               ##
##  Makefile for Building the LFSE library for Tinker/LFMM       ##
##                                                               ##
###################################################################
##
##  Invocation Options:
##
##   1. make all              Build all and make LSFE library
##   2. make clean            Delete objects and executables
##
##
###################################################################

# You want these to be consistent with those in the Makefile of Tinker

F77 = /usr/local/Cellar/gcc@7/7.5.0_4/bin/gfortran-7
LIBS =
F77FLAGS = -c
OPTFLAGS = -O
LIBFLAGS = -crusv
LINKFLAGS =

# This does not need to change unless you have moved this Makefile
TINKERDIR = ../source

###################################################################
##  No Change Needed Below This Line                             ##
###################################################################

OBJS = angsiz.o \
       caleul.o \
       clfse.o \
       cross.o \
       dcalc.o \
       ecalc.o \
       eulcnv.o \
       f01bcf.o \
       f02axf.o \
       f02ayf.o \
       fcalc.o \
       geom.o \
       hcross.o \
       iccglfse.o \
       matcal.o \
       p01aaf.o \
       power.o \
       print1.o \
       print2.o \
       rotate.o \
       sjan1d.o \
       sjmat.o \
       vcal.o \
       vdiag.o \
       wig3j.o \
       x02aaf.o \
       x02adf.o 

TNKLFMMLIB = libtinkerlfmm.a

TNKALLOBJ = $(wildcard $(TINKERDIR)/*.o)
TNKOBJ = $(filter-out %iccglfse.o, $(TNKALLOBJ))

%.o: %.f
	$(F77) $(F77FLAGS) $(OPTFLAGS) $< -o $@ 

all:    $(OBJS)
	ar ${LIBFLAGS} $(TNKLFMMLIB) $(TNKOBJ) $(OBJS)

clean:
	rm -f *.o *.a *.x

