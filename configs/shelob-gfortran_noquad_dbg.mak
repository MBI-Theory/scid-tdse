BUILD_ID :="Debug gfortran, built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*nq/    /' # Disable quad-math statements
F90 = gfortran-4.8.2 -I. -L/home/ps/lib64/ -Wl,-rpath,/home/ps/lib64/ \
      -m64 -mavx -g -Og -march=native -mtune=native -fopenmp \
      -std=f2003 -pedantic -Wall -Wno-unused-dummy-argument -Wno-unused-function \
      -Wno-maybe-uninitialized -Wcharacter-truncation -Wunderflow \
      -Wsurprising -fbacktrace -fcheck=all -fno-realloc-lhs \
      -ffpe-trap=invalid \
      -cpp -D__BUILD_ID__='$(BUILD_ID)' -ffree-line-length-none
F90L = $(F90) 
LAPACK = -llapack_gfortran -lopenblas 
LIBEXTRA = 

