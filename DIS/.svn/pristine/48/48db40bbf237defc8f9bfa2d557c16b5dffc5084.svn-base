#-*- Makefile -*-
###USER specific 
POWHEGPATH=$(PWD)/..
LHAPDF_CONFIG=$(shell which lhapdf-config)
FASTJET_CONFIG=$(shell which fastjet-config)
PYTHIA8_CONFIG=$(shell which pythia8-config)

## Choose compiler: gfortran,ifort
COMPILER=gfortran
## Choose Analysis: none, default
#ANALYSIS=sigmaR
#ANALYSIS=H1
#ANALYSIS=vbf
#ANALYSIS=caesar
#ANALYSIS=eic
#ANALYSIS=dummy
ANALYSIS=flux

ifeq ("$(COMPILER)","gfortran")	
F77= gfortran -fno-automatic 	-ffixed-line-length-none
CXX= g++  -std=c++11
## NOTE: -fbounds-check sometimes causes a compilation error
## due to non-lazy evaluation of boolean in gfortran.
# FFLAGS= -Wall -Wimplicit-interface
## For floating point exception trapping  uncomment the following 
#FPE=-ffpe-trap=invalid,zero,overflow,underflow 
## gfortran 4.4.1 optimized with -O3 yields erroneous results
## Use -O2 to be on the safe side
OPT=-O3
## For debugging uncomment the following
#DEBUG= -ggdb 
ifdef DEBUG
OPT=-O0
FPE=-ffpe-trap=invalid,zero,overflow
#,underflow
endif
endif

OBJ=obj-$(COMPILER)
OBJDIR:=$(OBJ)

STD=-lstdc++
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
STD=-lc++
endif

PWD=$(shell pwd)
WDNAME=$(shell basename $(PWD))
VPATH= ./:./POWHEGpatch:py8:$(POWHEGPATH):./amplitude-files/:$(OBJDIR)/

INCLUDE0=$(PWD)
INCLUDE1=$(PWD)/POWHEGpatch
INCLUDE2=$(POWHEGPATH)/include
INCLUDE3=amplitude-files
FF=$(F77) $(FFLAGS) $(FPE) $(OPT) $(DEBUG) -I$(INCLUDE0) -I$(INCLUDE1) -I$(INCLUDE2) -I$(INCLUDE3) 


INCLUDE =$(wildcard $(POWHEGPATH)/include/*.h *.h include/*.h amplitude-files/*.h)

PDFPACK=lhapdf6if.o lhapdf6ifcc.o
FJCXXFLAGS+= $(shell $(LHAPDF_CONFIG) --cxxflags)
LIBSLHAPDF= -Wl,-rpath,$(shell $(LHAPDF_CONFIG) --libdir)  -L$(shell $(LHAPDF_CONFIG) --libdir) -lLHAPDF $(STD)
LIBS+=$(LIBSLHAPDF)


ifeq ("$(ANALYSIS)","dummy")
	PWHGANAL=pwhg_analysis-dummy.o pwhg_bookhist-multi.o
else
	LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) $(STD)
	FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
	PWHGANAL=fastjetfortran.o pwhg_bookhist-multi.o 
	ifeq ("$(ANALYSIS)","H1")
		PWHGANAL+=pwhg_analysis_H1_2006.o event_shapes.o
	endif	
	ifeq ("$(ANALYSIS)","vbf")
		PWHGANAL+=pwhg_analysis_vbf.o 
	endif
	ifeq ("$(ANALYSIS)","sigmaR")
		PWHGANAL+=pwhg_analysis_inclusive_sigma.o 
	endif
	ifeq ("$(ANALYSIS)","caesar")
		PWHGANAL+=pwhg_analysis_caesar.o event_shapes.o
	endif
	ifeq ("$(ANALYSIS)","eic")
		PWHGANAL+=pwhg_analysis_eic.o 
	endif
	ifeq ("$(ANALYSIS)","flux")
		PWHGANAL+=pwhg_analysis_flux.o 
	endif
endif

PWHGANAL+=multi_plot.o 

%.o: %.f $(INCLUDE) | $(OBJDIR)
	$(FF) -c -o $(OBJDIR)/$@ $<

%.o: %.f90 $(INCLUDE) | $(OBJDIR)
	$(FF)  -ffree-line-length-512 -c -o $(OBJ)/$@ $<

%.o: %.c | $(OBJDIR)
	$(CC) $(DEBUG) -c -o $(OBJ)/$@ $^ 

%.o: %.cc | $(OBJDIR)
	$(CXX) $(DEBUG) -c -o $(OBJ)/$@ $^ $(FJCXXFLAGS)  $(LIBSFASTJET)


LIBSAMP=-Lamplitude-files -lamplitude 
LIBS+=-lz 
USER=init_couplings.o init_processes.o Born_phsp.o Born.o virtual.o	\
     real.o lepton_flux.o $(PWHGANAL)

AMP=braket_vbf.o help_hh.o qq_ee.o qqj_ee.o

PWHG=pwhg_main.o pwhg_init.o bbinit.o btilde.o lhefwrite.o		\
	LesHouches.o LesHouchesreg.o gen_Born_phsp.o find_regions.o	\
	test_Sudakov.o pt2maxreg.o sigborn.o gen_real_phsp.o maxrat.o	\
	gen_index.o gen_radiation.o Bornzerodamp.o sigremnants.o	\
	random.o boostrot.o cernroutines.o	\
	init_phys.o powheginput.o pdfcalls.o sigreal.o sigcollremn.o	\
	pwhg_analysis_driver.o checkmomzero.o		\
	setstrongcoupl.o integrator.o newunit.o mwarn.o sigsoftvirt.o	\
	sigcollsoft.o sigvirtual.o reshufflemoms.o ubprojections-new.o      \
	setlocalscales.o mint_upb.o opencount.o          \
	validflav.o $(PDFPACK) $(USER) $(FPEOBJ) lhefread.o pwhg_io_interface.o \
	rwl_weightlists.o rwl_setup_param_weights.o rwl_setup_param_weights_user.o \
	mintwrapper.o fullrwgt.o locks.o cache_similar.o sigequiv_hook.o fill_res_histories.o \
	resweights.o sigregular.o genericphsp.o boostrot4.o PhaseSpaceUtils.o\
	convert_coup.o utils.o bra_ket_subroutines.o elliptc_fns.o 


# target to generate LHEF output
pwhg_main: amplitude-files/libamplitude.a $(PWHG)
	$(FF) $(patsubst %,$(OBJ)/%,$(PWHG)) $(LIBS) $(LIBSAMP) $(LIBSFASTJET) $(STATIC) -o $@

amplitude-files/libamplitude.a: amplitude-files/*f 
	cd amplitude-files ; make


# target to generate object directory if it does not exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

LHEF=lhef_analysis.o boostrot.o random.o cernroutines.o	\
     opencount.o powheginput.o $(PWHGANAL)	\
     lhefread.o pwhg_io_interface.o rwl_weightlists.o newunit.o pwhg_analysis_driver.o locks.o $(FPEOBJ)

# target to analyze LHEF output
lhef_analysis: $(LHEF) 
	$(FF) $(patsubst %,$(OBJ)/%,$(LHEF)) $(LIBS) $(LIBSFASTJET) $(STATIC) -o $@


# PYTHIA 8
FJCXXFLAGS+=$(shell $(PYTHIA8_CONFIG) --cxxflags) $(STD)
LIBPYTHIA8=$(shell $(PYTHIA8_CONFIG) --ldflags) 
LIBHEPMC=-lHepMC

# target to read event file, shower events with PYTHIA8 + analysis
PYTHIA8=main-PYTHIA8.o pythia83F77.o powheginput_cpp.o	      \
	powheginput.o pwhg_io_interface.o rwl_weightlists.o	      \
	$(PWHGANAL) opencount.o lhefread.o newunit.o pdfdummies.o     \
	reshufflemoms.o random.o cernroutines.o bra_ket_subroutines.o \
	Check_LesHouches.o locks.o utils.o boostrot.o $(FPEOBJ) $(LIBZDUMMY)


main-PYTHIA8-lhef: $(PYTHIA8) 
	$(FF) $(LIBPYTHIA8) $(patsubst %,$(OBJ)/%,$(PYTHIA8)) $(LIBSFASTJET) $(LIBPYTHIA8) $(LIBHEPMC) $(RIVETLIBS) $(STATIC) $(LIBLHAPDF) $(LIBS) -o $@



# target to cleanup
.PHONY: clean

clean:
	rm -f $(OBJ)/*.o pwhg_main lhef_analysis	\
	main-PYTHIA8-lhef amplitude-files/*.o amplitude-files/libamplitude.a

