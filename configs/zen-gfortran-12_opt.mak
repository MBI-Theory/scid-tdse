BUILD_ID :="Optimized gfortran-12 (ZEN/EPYC), built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
F90 = gfortran-12 -I. \
      -flto -O3 -fprotect-parens -march=native -mtune=native -fopenmp \
      -ffast-math -fcx-fortran-rules -mrecip \
      -floop-block \
      -fno-realloc-lhs -fbacktrace -g \
      -static-libgfortran \
      -cpp -D__BUILD_ID__='$(BUILD_ID)' -ffree-line-length-none
#     -fexternal-blas -fblas-matmul-limit=50 # Causes memory errors
#     -floop-block # Requires Graphite isl, not configured on many distributions
#     -ffast-math -fcx-fortran-rules -mrecip 
F90L = $(F90) 
LAPACK = -lopenblas_openmp # -lquadlapack_gfortran
LIBEXTRA = -B/usr/share/libhugetlbfs -lhugetlbfs -Wl,--hugetlbfs-align
