BUILD_ID :="Debug gfortran, built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
F90 = gfortran -I. \
      -m64 -mavx -g -Og -march=native -mtune=native -fopenmp \
      -std=gnu -pedantic -Wall -Wno-unused-dummy-argument -Wno-unused-function \
      -Wno-maybe-uninitialized -Wcharacter-truncation -Wunderflow \
      -Wsurprising -fbacktrace -fcheck=all -fno-realloc-lhs \
      -ffpe-trap=invalid \
      -cpp -D__BUILD_ID__='$(BUILD_ID)' -ffree-line-length-none
F90L = $(F90) 
# NB liblapack coming with OpenSuse is broken for OpenMP!
LAPACK = /home/ps/lib64/liblapack_gfortran.a -lopenblas_openmp # -lquadlapack_gfortran
LIBEXTRA = 

