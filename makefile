.PHONY: clean 

OMP_STACKSIZE        := 8192M
export OMP_STACKSIZE

OMP_NUM_THREADS      := $(shell nproc)
export OMP_NUM_THREADS

# OMP_NESTED           := true
# export OMP_NESTED


# -fdefault-real-8 should be added to FPREC to promote doubles to REAL*16
FCOMP     = gfortran-7
FPREC     = -freal-8-real-10
FOPTS     = -w -O2 -march=native -fopenmp  $(FPREC)
LDOPT     = 
PYTHON    = python


FC     = $(FCOMP)
FFLAGS = $(FOPTS)
export FPREC

CC        = gcc
COPTS     = -I./include -O3

# Set the list of programs to compile

PROGRAMS               = test_chebyshev test_riccati test_pbessel test_helmrad            \
                         test_helmop expr_pbessel expr_pbessel2 expr_pbessel3             \
                         expr_pbessel4 expr_bump expr_volcano expr_discont                

##########################################################################################



EXPR_DISCONT_FILES     = utils.o                                                         \
                         chebyshev.o                                                     \
                         riccati.o                                                       \
                         hank103.o                                                       \
                         pbessel.o                                                       \
                         helmrad.o                                                       \
                         helmop.o                                                        \
                         libdfftpack.a

EXPR_VOLCANO_FILES     = utils.o                                                         \
                         chebyshev.o                                                     \
                         riccati.o                                                       \
                         hank103.o                                                       \
                         pbessel.o                                                       \
                         helmrad.o                                                       \
                         helmop.o                                                        \
                         libdfftpack.a

EXPR_BUMP_FILES        = utils.o                                                         \
                         chebyshev.o                                                     \
                         hank103.o                                                       \
                         riccati.o                                                       \
                         pbessel.o                                                       \
                         helmrad.o                                                       \
                         helmop.o                                                        \
                         libdfftpack.a

EXPR_PBESSEL_FILES     = utils.o                                                         \
                         chebyshev.o                                                     \
                         riccati.o                                                       \
                         pbessel.o                                                       \

HELMOP_FILES           = utils.o                                                         \
                         chebyshev.o                                                     \
                         helmop.o                                                        \
                         libdfftpack.a

HELMRAD_FILES          = utils.o                                                         \
                         chebyshev.o                                                     \
                         hank103.o                                                       \
                         riccati.o                                                       \
                         pbessel.o                                                       \
                         helmrad.o                                                       \
                         libdfftpack.a

PBESSEL_FILES          = utils.o                                                         \
                         chebyshev.o                                                     \
                         riccati.o                                                       \
                         pbessel.o                                                       \

RICCATI_FILES          = utils.o                                                         \
                         chebyshev.o                                                     \
                         riccati.o

CHEBYSHEV_FILES        = utils.o                                                         \
                         chebyshev.o

##########################################################################################

all	                : clean $(PROGRAMS) 


libdfftpack.a           :
	cd dfftpack && $(MAKE) clean && $(MAKE) && cp libdfftpack.a ..

mtest.o                 : $(MTEST_FILES) mtest.f90
mtest                   : $(MTEST_FILES) mtest.o

expr_discont.o          : $(EXPR_DISCONT_FILES) expr_discont.f90
expr_discont            : $(EXPR_DISCONT_FILES) expr_discont.o

expr_volcano.o          : $(EXPR_VOLCANO_FILES) expr_volcano.f90
expr_volcano            : $(EXPR_VOLCANO_FILES) expr_volcano.o

expr_bump               : $(EXPR_BUMP_FILES) expr_bump.o
expr_bump.o             : $(EXPR_BUMP_FILES) expr_bump.f90

expr_pbessel4.o         : $(EXPR_PBESSEL_FILES) expr_pbessel4.f90
expr_pbessel4           : $(EXPR_PBESSEL_FILES) expr_pbessel4.o

expr_pbessel3.o         : $(EXPR_PBESSEL_FILES) expr_pbessel3.f90
expr_pbessel3           : $(EXPR_PBESSEL_FILES) expr_pbessel3.o

expr_pbessel2.o         : $(EXPR_PBESSEL_FILES) expr_pbessel2.f90
expr_pbessel2           : $(EXPR_PBESSEL_FILES) expr_pbessel2.o

expr_pbessel.o          : $(EXPR_PBESSEL_FILES) expr_pbessel.f90
expr_pbessel            : $(EXPR_PBESSEL_FILES) expr_pbessel.o

test_helmop.o           : $(HELMOP_FILES) test_helmop.f90
test_helmop             : $(HELMOP_FILES) test_helmop.o

test_helmrad.o          : $(HELMRAD_FILES) test_helmrad.f90
test_helmrad            : $(HELMRAD_FILES) test_helmrad.o

test_pbessel.o          : $(PBESSEL_FILES) test_pbessel.f90
test_pbessel            : $(PBESSEL_FILES) test_pbessel.o

test_odesolve.o         : $(ODESOLVE_FILES) test_odesolve.f90
test_odesolve           : $(ODESOLVE_FILES) test_odesolve.o

test_riccati.o          : $(RICCATI_FILES) test_riccati.f90
test_riccati            : $(RICCATI_FILES) test_riccati.o

test_chebyshev.o        : $(CHEBYSHEV_FILES) test_chebyshev.f90
test_chebyshev          : $(CHEBYSHEV_FILES) test_chebyshev.o



# Setup the general compilation rules

%		: %.o
	$(FCOMP) $(FOPTS) -o $@ $^ $(LDOPT)
	@echo  
	@echo 
	@echo "---------[ $@     ]--------------------------------------------"
	@echo 
	@./$@
	@echo 
	@echo "--------------------------------------------------------------------------"
	@echo 

%.o		: %.f90
	$(FCOMP) -c $(FOPTS)  $<

%.o		: %.f90
	$(FCOMP) -c $(FOPTS)  $<

%.o		: %.f
	$(FCOMP) -c $(FOPTS)  $<

%.o		: %.cpp
	$(CPPCOMP) -c $(CPPOPTS)  $<

%.o		: %.cpp
	$(CPPCOMP) -c $(CPPOPTS)  $<

%.o		: %.c
	$(CC) -c $(COPTS)  $<

.PHONY: clean

clean:
	rm -f *.o *.mod *~ fort.* a.out *.a
	rm -f $(PROGRAMS)
	rm -f *.py
	rm -f *.pdf
	rm -f *.tex
	cd dfftpack; make clean; cd ..


