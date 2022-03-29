# start of the makefile
# defining variables
f90comp = gfortran

# unknown starts

FFLAGS =

#PLIBS = -L/projects/intel/composerxe-2011.1.107/composerxe-2011.5.220/mkl/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -L/home/hp5/Sources/gsl/lib -lgsl
#PLIBS = -L/opt/intel/composer_xe_2013.1.117/mkl/lib/intel64 -lmkl_lapack95_lp64  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -L/usr/lib64 -lpthread
#PLIBS = -L/usr/lib64/atlas -llapack -lf77blas -lcblas -latlas
#PLIBS = -L/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -L/usr/lib64 -lpthread
#LIBS = $(PLIBS)
MKLROOT=/gpfsnyu/packages/intelPSX/2019/compilers_and_libraries_2019.0.117/linux/mkl

PLIBS = -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
LIBS = $(PLIBS)

# unknown ends

# DMFT_chi_multiple: DMFT_chi_multiple.o do_chi_multiple.o
# $(f90comp) $(FFLAGS) -o DMFT_chi_multiple DMFT_chi_multiple.o do_chi_multiple.o $(LIBS)

# makefile
power: power.o
	$(f90comp) -llapack -mcmodel=large power.o -o power

power.o: power.f
	$(f90comp) -c power.f

%.o: %.f
	$(f90comp) -c $<

# cleaning everything
clean:
	rm -rf power.o power
# end of the makefile

