BUILD_ID :="Optimized gfortran, built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
# scl run gcc-toolset-14 bash
# -m64 -mavx -O3 -march=core-avx-i -mtune=core-avx-i -fopenmp
# Zen v2 binaries also work on AVX2-capable Intel CPUs
F90 = /opt/rh/gcc-toolset-14/root/usr/bin/gfortran -I. -I${MKLROOT}/include \
      -std=gnu -Wall -Wcharacter-truncation -Wextra -Wno-unused-dummy-argument -Wno-compare-reals \
         -Wimplicit-interface -Wimplicit-procedure -Wintrinsics-std \
         -Wrealloc-lhs-all \
      -flto -O3 -fprotect-parens -march=znver2 -mtune=znver2 -fopenmp \
      -ffast-math -fcx-fortran-rules -freciprocal-math \
      -floop-block \
      -fno-realloc-lhs -fbacktrace -g \
      -static-libgfortran -static-libquadmath \
      -cpp -D__BUILD_ID__='$(BUILD_ID)' -ffree-line-length-none
#     -ffast-math -fcx-fortran-rules -mrecip
F90L = $(F90) 
#LAPACK = -L/home/ps/lib64 -llapack_gfortran -lblas 
LAPACK =  -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_gf_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group -lpthread -lm -ldl
LIBEXTRA = -B/usr/share/libhugetlbfs -lhugetlbfs -Wl,--hugetlbfs-align
