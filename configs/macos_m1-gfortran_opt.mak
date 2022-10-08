BUILD_ID :="Optimized gfortran, built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*nq/    /' # Disable quad-math statements
# -Ofast is known to produce code which crashes on M1, at least with gfortran 12.2
F90 = gfortran -I. \
      -flto -O1 -fprotect-parens -march=native -mtune=native -fopenmp \
      -ffast-math -fcx-fortran-rules \
      -floop-block \
      -fno-realloc-lhs -fbacktrace -g \
      -cpp -D__BUILD_ID__='$(BUILD_ID)' -ffree-line-length-none
F90L = $(F90) 
#
#  Both openblas and Accelerate are known to work, with openblas being somewhat faster
#
#LAPACK = -L/opt/homebrew/opt/openblas/lib -lopenblas
LAPACK = -framework Accelerate
LIBEXTRA = 
