BUILD_ID :="Optimized ifort, built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
# -xAVX is broken in 16.0 and later versions of 15.0. As usual.
# The latest version known to work is 15.0.1.133
F90 = ifort \
            -mmic -I$(MKLROOT)/include \
            -fast \
            -openmp -align all -align array256byte -pad \
            -warn -assume buffered_io \
            -complex_limited_range -ftz -heap-arrays 32 \
            -qopt-prefetch=0 -opt-matmul -opt-subscript-in-range -opt-dynamic-align \
            -qopt-mem-layout-trans=3 -qopt-multi-version-aggressive \
            -debug extended -traceback \
            -cpp -D__BUILD_ID__='$(BUILD_ID)'
F90L = $(F90) 
LAPACK = # -L/home/ps/lib64 -lquadlapack_intel
# -lhugetlbfs -Wl,-z,common-page-size=2097152 -Wl,-z,max-page-size=2097152
# LIBEXTRA =  -Wl,--start-group $(MKLROOT)/lib/mic/libmkl_intel_lp64.a $(MKLROOT)/lib/mic/libmkl_core.a $(MKLROOT)/lib/mic/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl
LIBEXTRA =  -Wl,--start-group $(MKLROOT)/lib/mic/libmkl_intel_lp64.a $(MKLROOT)/lib/mic/libmkl_core.a $(MKLROOT)/lib/mic/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm -ldl
