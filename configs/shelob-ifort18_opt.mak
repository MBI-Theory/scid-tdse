BUILD_ID :="Optimized ifort, built on $(shell hostname) at $(shell date)"
ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
# -xAVX is broken in 16.0 and later versions of 15.0. As usual.
# The latest version known to work is 15.0.1.133
F90 = ifort \
            -qopenmp -align all -align array256byte -pad \
            -warn -assume buffered_io \
            -O3 -ipo8 -no-prec-div -xAVX -axAVX2,MIC-AVX512 -complex_limited_range -fp-model fast=1 -ftz -assume protect_parens -heap-arrays 32 \
            -qopt-prefetch=0 -qopt-matmul -qopt-subscript-in-range -qopt-dynamic-align \
            -qopt-mem-layout-trans=3 -qopt-multi-version-aggressive \
            -debug full -debug extended -traceback -qopt-report=5 \
            -cpp -D__BUILD_ID__='$(BUILD_ID)'
#  Conservative: -complex_limited_range -fp-model fast=1 -ftz -assume protect_parens
#  Aggressive:   -complex_limited_range -fp-model fast=2 -ftz
# -fp-model precise -assume protect_parens # is necessary to make ifort fully compliant with the Fortran semantics
# -complex_limited_range -fp-model strict -ftz # produces code similar to gfortran in speed, but still less accurate
# -complex_limited_range -fp-model fast=1 -ftz # produces faster code, but dramatically decreases the accuracy
# -complex_limited_range -fp-model fast=2 -ftz # produces still slightly faster code, but dramatically decreases the accuracy
F90L = $(F90) 
LAPACK = # -L/home/ps/lib64 -lquadlapack_intel
LIBEXTRA = -lhugetlbfs -Wl,-z,common-page-size=2097152 -Wl,-z,max-page-size=2097152 \
           -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm \
           -Wl,-rpath,$(MKLROOT)/lib/intel64/ -Wl,-rpath,/opt/intel/compilers_and_libraries_2018.0.128/linux/compiler/lib/intel64/
