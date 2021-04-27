BUILD_ID :="Optimized ifort, built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
# -xAVX is broken in 16.0 and later versions of 15.0. As usual.
# The latest version known to work is 15.0.1.133
F90 = ifort \
            -openmp -align all -align array256byte -pad \
            -warn -assume buffered_io \
            -O3 -ipo -no-prec-div -xCORE-AVX2 -complex_limited_range -fp-model fast=2 -ftz -heap-arrays 32 \
            -qopt-prefetch=0 -opt-matmul -opt-subscript-in-range -opt-dynamic-align \
            -qopt-mem-layout-trans=3 -qopt-multi-version-aggressive \
            -debug full -debug extended -traceback -opt-report=5 \
            -cpp -D__BUILD_ID__='$(BUILD_ID)'
F90L = $(F90) 
LAPACK = # -L/home/ps/lib64 -lquadlapack_intel
LIBEXTRA = -lhugetlbfs -Wl,-z,common-page-size=2097152 -Wl,-z,max-page-size=2097152 \
           -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm \
           -Wl,-rpath,$(MKLROOT)/lib/intel64/ -Wl,-rpath,/opt/intel/composer_xe_2015.1.133/compiler/lib/intel64/
