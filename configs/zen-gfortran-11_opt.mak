BUILD_ID :="Optimized gfortran-11 (ZEN/EPYC), built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
# WARNING: requesting -mtune=znver2 leads to incorrect code
# WARNING: requesting -fexternal-blas leads to unreasonable stack requirements
F90 = gfortran-11 -I. \
      -O3 -fprotect-parens -march=native -mtune=znver1 -fopenmp \
      -ffast-math -fcx-fortran-rules -mrecip \
      -fno-realloc-lhs -fbacktrace -g \
      -static-libgfortran \
      -cpp -D__BUILD_ID__='$(BUILD_ID)' -ffree-line-length-none
#     -fexternal-blas -fblas-matmul-limit=50 
#     -floop-block # Requires Graphite isl, not configured on many distributions
#     -ffast-math -fcx-fortran-rules -mrecip 
F90L = $(F90) 
LAPACK = -lopenblas_openmp # -lquadlapack_gfortran
LIBEXTRA = -B/usr/share/libhugetlbfs -lhugetlbfs -Wl,--hugetlbfs-align
