FC = gfortran
POWHEGPATH=$(PWD)/../..
STD=-lstdc++
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
STD=-lc++
endif

INCLUDE1=$(PWD)
INCLUDE2=$(PWD)/..
INCLUDE3=$(PWD)/../POWHEGpatch/
INCLUDE4=$(POWHEGPATH)/include

ifeq ("$(FC)","gfortran")	
FFLAGS = -g -fno-automatic -ffixed-line-length-none -O2 
endif
ifeq ("$(FC)","ifort")	
FFLAGS = -save -extend_source -O3  
endif
FFLAGS+= -I$(INCLUDE1) -I$(INCLUDE2) -I$(INCLUDE3)  -I$(INCLUDE4) 	

# ----------------------------------------------------------------------------

# object files for qq->qqHH

OBJECTS = ./qq_ee.o \
          ./qqj_ee.o \
          ./help_hh.o \
	  ./braket_vbf.o  

TARGETS = libamplitude.a


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# rules for vbfnlo_lib:
#
 libamplitudes.a:  $(OBJECTS)
	$(AR) r libamplitude.a $(OBJECTS) 

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

# for g77 or gfortran: 
.f.o:
	$(FC) -c -o $*.o  $(FFLAGS) $*.f

# ----------------------------------------------------------------------------

clean:
	rm -f *.o core *~ *.x libamplitude.a



