BUILD_ID :="Optimized gfortran, built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*nq/    /' # Disable quad-math statements
F90 = gfortran-4.8.2 -I. -L/home/ps/lib64/ -Wl,-rpath,/home/ps/lib64/ \
      -m64 -mavx -O3 -march=native -mtune=native -fopenmp \
      -floop-block \
      -ffast-math -fcx-fortran-rules -mrecip \
      -fexternal-blas -fblas-matmul-limit=50 \
      -fno-realloc-lhs -fbacktrace -g \
      -cpp -D__BUILD_ID__='$(BUILD_ID)' -ffree-line-length-none
F90L = $(F90) 
LAPACK = -llapack_gfortran -lopenblas
LIBEXTRA = -B/usr/share/libhugetlbfs -lhugetlbfs -Wl,--hugetlbfs-align

