BUILD_ID :="Optimized gfortran, built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*nq/    /' # Enable quad-math statements
# WARNING: gfortran 11.2.0 appears to have many bugs in OpenMP code generation.
# WARNING: Please make sure that the ENTIRE test suite runs correctly before using in production.
F90 = gfortran -I. \
      -O3 -march=native -mtune=native -fopenmp \
      -floop-block \
      -ffast-math -fcx-fortran-rules \
      -fno-realloc-lhs -fbacktrace -g \
      -cpp -D__BUILD_ID__='$(BUILD_ID)' -ffree-line-length-none
F90L = $(F90) 
LAPACK = -L/opt/homebrew/opt/openblas/lib -lopenblas
LIBEXTRA = 
