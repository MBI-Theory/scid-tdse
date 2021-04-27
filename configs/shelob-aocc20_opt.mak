#
#  Don't forget: "source /opt/AMD/aocc-compiler-2.3.0/setenv_AOCC.sh"
#  It may be necessary to edit the defintion of isnan in accuracy.f90
#
#  WARNING: 
#
#  As of version 2.3.0, AOCC's flang generates amateur-quality code
#  code, with gfortran and ifort outperforming it for all test cases,
#  usually by a large margin (2x or so). Note that AOCC miscompiles
#  timer.f90, so that the CPU times cannot be reported. The real-time
#  clock provided by the system_clock() intrinsic stops counting after
#  2147.5 seconds.
#
#  The use of this "compiler" is STRONGLY DISCOURAGED.
#
BUILD_ID :="Optimized AOCC 2.3.0, built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*nq/    /' # Disable quad-math statements; AOCC apparently does not support them
# Flags from peak SPEC submission by Supermicro
F90 = flang \
           -fopenmp \
           -O3 -flto -march=native -funroll-loops -Mrecursive -mllvm -vector-library=LIBMVEC \
           -Kieee -fno-finite-math-only \
           -cpp -D__BUILD_ID__='$(BUILD_ID)'
#
timer.o: timer.f90
	$(ACT) $(ACT2) $< >preprocess/$<
	flang -fopenmp -O0 -Mrecursive -cpp -D__BUILD_ID__='$(BUILD_ID)' -c -o timer.o preprocess/$<
#
F90L = $(F90) \
           -Wl,-mllvm -Wl,-function-specialize \
           -Wl,-mllvm -Wl,-region-vectorize \
           -Wl,-mllvm -Wl,-vector-library=LIBMVEC \
           -Wl,-mllvm -Wl,-reduce-array-computations=3
LAPACK = /opt/AMD/aocl/aocl-linux-aocc-2.2.0/lib/libflame.a \
         /opt/AMD/aocl/aocl-linux-aocc-2.2.0/lib/libblis.a \
         /opt/AMD/aocc-compiler-2.3.0/lib/libamdlibm.a
#LAPACK = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_gf_lp64.a \
#                           $(MKLROOT)/lib/intel64/libmkl_core.a \
#                           $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -ldl
LIBEXTRA = -lhugetlbfs -Wl,-z,max-page-size=2097152 \
           -lmvec -ljemalloc -lm 
