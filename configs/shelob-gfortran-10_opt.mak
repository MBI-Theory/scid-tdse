BUILD_ID :="Optimized gfortran-10, built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
# Add m64 -mavx -mavx2 etc non-native instruction sets are wanted.
F90 = gfortran-10 -I. \
      -O3 -flto -fprotect-parens -march=native -mtune=native -fopenmp \
      -ffast-math -fcx-fortran-rules -mrecip \
      -fno-realloc-lhs -fbacktrace -g \
      -static-libgfortran \
      -cpp -D__BUILD_ID__='$(BUILD_ID)' -ffree-line-length-none
#
#     -fexternal-blas -fblas-matmul-limit=50 # fucked up in gfortran 10
#     -floop-block # Requires Graphite isl, not configured on many distributions
F90L = $(F90) 
#LAPACK = -lopenblas_openmp
LAPACK = /home/ps/lib64/liblapack_gfortran.a -lopenblas_openmp # -lquadlapack_gfortran
LIBEXTRA = -B/usr/share/libhugetlbfs -lhugetlbfs -Wl,--hugetlbfs-align
