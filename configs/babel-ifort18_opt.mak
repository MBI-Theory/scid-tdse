BUILD_ID :="Optimized ifort18, built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
# 
F90 = /opt/intel/oneapi/compiler/2021.1.1/linux/bin/intel64/ifort \
            -qopenmp -align all -align array256byte -pad \
            -warn -assume buffered_io \
            -O3 -ipo -no-prec-div -xAVX -axAVX2,MIC-AVX512 -complex_limited_range -fp-model fast=1 -ftz -assume protect_parens -heap-arrays 32 \
            -qopt-prefetch=0 -qopt-matmul -qopt-subscript-in-range -qopt-dynamic-align \
            -qopt-mem-layout-trans=3 -qopt-multi-version-aggressive \
            -debug full -debug extended -traceback -qopt-report=5 \
            -cpp -D__BUILD_ID__='$(BUILD_ID)'
F90L = $(F90) 
LAPACK = 
LIBEXTRA = -lhugetlbfs -Wl,-z,common-page-size=2097152 -Wl,-z,max-page-size=2097152 \
           -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm \
           -Wl,-rpath,$(MKLROOT)/lib/intel64/ -Wl,-rpath,/opt/intel/oneapi/compiler/2021.1.1/linux/compiler/lib/intel64_lin/
