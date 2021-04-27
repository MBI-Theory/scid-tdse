BUILD_ID :="Optimized gfortran-7 (ZEN/EPYC), built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
F90 = gfortran-7 -I. \
      -m64 -mavx2 -O3 -fprotect-parens -march=znver1 -mtune=znver1 -fopenmp \
      -ffast-math -fcx-fortran-rules -mrecip \
      -fexternal-blas -fblas-matmul-limit=50 \
      -fno-realloc-lhs -fbacktrace -g \
      -static-libgfortran \
      -cpp -D__BUILD_ID__='$(BUILD_ID)' -ffree-line-length-none
#     -floop-block # Requires Graphite isl, not configured on many distributions
#     -ffast-math -fcx-fortran-rules -mrecip 
F90L = $(F90) 
LAPACK = /home/ps/lib64/liblapack_gfortran.a -lopenblas # -lquadlapack_gfortran
LIBEXTRA = -B/usr/share/libhugetlbfs -lhugetlbfs -Wl,--hugetlbfs-align
