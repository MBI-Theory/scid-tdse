BUILD_ID :="Optimized Intel oneAPI (8-byte integers), built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
#ACT2 = -e 's/^!\*mp/    /' # Enable MPI statements
# -axMIC-AVX512 is not supported by recent oneapi
F90 = ifort \
            -i8 -qopenmp -qmkl=sequential -align all -align array256byte -pad \
            -warn -assume buffered_io \
            -O3 -ipo8 -no-prec-div -xAVX -axAVX2 -complex_limited_range -fp-model fast=1 -ftz -assume protect_parens -heap-arrays 32 \
            -qopt-prefetch=0 -qopt-matmul -qopt-subscript-in-range -qopt-dynamic-align \
            -qopt-mem-layout-trans=3 -qopt-multi-version-aggressive \
            -debug full -debug extended -traceback -qopt-report=5 \
            -cpp -D__BUILD_ID__='$(BUILD_ID)'
F90L = $(F90) 
LAPACK = # Lapack is in Intel MKL, 
LIBEXTRA = \
           -static-intel -qopenmp-link=static \
           -lhugetlbfs -Wl,-z,common-page-size=2097152 -Wl,-z,max-page-size=2097152 \
           -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm 
