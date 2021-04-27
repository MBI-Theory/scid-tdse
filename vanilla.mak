#
#  Default build configuration file. 
#
#  See README.txt for further instructions.
#
BUILD_ID :="Plain vanilla, built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*nq/    /' # Disable quad-math statements
#  Fortran compiler. These optimization flags are most likely grossly inadequate;
#  see examples in the configs/ subdirectory for suggestions.
F90 = gfortran -O -cpp -I. -D__BUILD_ID__='$(BUILD_ID)'
#  Fortran linker
F90L = $(F90) 
#  BLAS and LAPACK libraries
#  If you enable OpenMP parallel execution, the libraries you supply MUST support
#  multi-threaded execution as well!
LAPACK = -llapack -lblas
#  Any additional support libraries
LIBEXTRA = 
