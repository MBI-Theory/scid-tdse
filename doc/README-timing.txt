Last updated: 2021 Feb 17
-------------

Executing the test cases in ./examples (which is strongly recommended
in any case, as dodgy libraries and compiler bugs are, sadly, very
common!) will also produce timing data. As a side effect, it also
produces timing information (use timing-summary.awk to extract it from
the output of run-all.sh; see however note [6] below). 

As a reference point, some timing data is collected below.

Run times [in seconds] for the test cases in ./examples. 

 Test case                             gf7[1,2,3]  gf10[1,2,4]  if21[1,2,5]  aoc23[1,2,6]
 ------------                          ----------  -----------  -----------  ------------
 hydrogen_1S_2P0_uniform.inp                0.8          0.8          1.6          2.4
 hydrogen_1S_2P0_uniform_restart.inp        0.8          0.7          1.0          2.1
 hydrogen_1S_2P0.inp                       25.1         22.0         23.8         48.8
 hydrogen_2P0_ion.inp                      35.9         32.0         32.0         68.6
 hydrogen_2P0_ion_restart.inp              24.9         25.1         18.2         52.2
 helium_GJG75.inp                          48.7         43.8         43.2         95.6
 helium_triplet_spline_linear.inp          27.1         24.7         16.8         41.3
 helium_1S_adiabatic.inp                  107.9         92.4        115.3        196.5
 argon_3P1_cooper.inp                     107.3        102.1         97.0        199.8
 argon_3P1_offcooper.inp                  118.6        114.7        104.1        212.8
 argon_3P1m_circ_l.inp                    147.5        151.5        112.0        297.5
 argon_3P1m_circ_r.inp                    312.2        274.2        244.9        679.5
 hydrogen_1S_hhg_linear.inp               211.1        185.5        230.6        305.9
 hydrogen_2P0_sfi_tsurf.inp               603.2        551.1        506.8        877.9
 argon_3P1m_ell_ckpt_mpi.inp              321.8        257.0        285.1        501.7
 argon_3P1m_ell_rstrt_mpi.inp             109.3         86.7         96.3        168.9
 hydrogen_1S_hhg_spline.inp               253.5        226.8        213.3        482.3
 hydrogen_2P0_sfi.inp                    5061.2       4352.8       4663.3       7388.
 hydrogen_1S_hhg_elliptical.inp          2139.2       1725.8       2006.0       3496.
 hydrogen_1S_hhg_circular.inp            4100.2       3135.2       3847.3       6267.
 hydrogen_1S_hhg_ell_refsol.inp          4914.6       4216.0       4413.5       7112.
 ------------                           -------      -------      -------      -------
 Total                                  18670.9      15620.9      17072.1      28497.

[1] Source versions used in the timing test:

    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $
    $Id: README-timing.txt,v 1.1 2021/02/17 16:58:44 ps Exp $

[2] AMD Ryzen Threadripper 3960X 24-Core Processor, 4-channel DDR3200 (non-ECC) RAM

[3] gfortran 7.5.0, with the following flags: 
      -I. -m64 -mavx -O3 -fprotect-parens -march=native -mtune=native -fopenmp \
      -ffast-math -fcx-fortran-rules -mrecip -fexternal-blas -fblas-matmul-limit=50 \
      -fno-realloc-lhs -fbacktrace -g -static-libgfortran \
      -cpp -D__BUILD_ID__='$(BUILD_ID)' -ffree-line-length-none \
      /home/ps/lib64/liblapack_gfortran.a -lopenblas_openmp \
      -B/usr/share/libhugetlbfs -lhugetlbfs -Wl,--hugetlbfs-align

    Minor failures in hydrogen_2P0_sfi.inp and hydrogen_1S_hhg_circular.inp

[4] gfortran 10.2.1, with the following flags:
      -I. -O3 -flto -fprotect-parens -march=native -mtune=native -fopenmp \
      -ffast-math -fcx-fortran-rules -mrecip -fno-realloc-lhs -fbacktrace -g \
      -static-libgfortran -cpp -D__BUILD_ID__='$(BUILD_ID)' -ffree-line-length-none \
      /home/ps/lib64/liblapack_gfortran.a -lopenblas_openmp \
      -B/usr/share/libhugetlbfs -lhugetlbfs -Wl,--hugetlbfs-align

    Minor failures in hydrogen_2P0_sfi.inp, hydrogen_1S_hhg_circular.inp, and
    hydrogen_1S_hhg_spline.inp

[5] Intel Fortran 2021.1 Build 20201112_000000, with the following flags:
      -qopenmp -align all -align array256byte -pad -warn -assume buffered_io \
      -O3 -ipo8 -no-prec-div -xAVX -axAVX2,MIC-AVX512 -complex_limited_range \
      -fp-model fast=1 -ftz -assume protect_parens -heap-arrays 32 \
      -qopt-prefetch=0 -qopt-matmul -qopt-subscript-in-range -qopt-dynamic-align \
      -qopt-mem-layout-trans=3 -qopt-multi-version-aggressive \
      -debug full -debug extended -traceback -qopt-report=5 \
      -cpp -D__BUILD_ID__='$(BUILD_ID)' \
      -lhugetlbfs -Wl,-z,common-page-size=2097152 -Wl,-z,max-page-size=2097152 \
      -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
      $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a \
      -Wl,--end-group -lpthread -lm \
      -Wl,-rpath,$(MKLROOT)/lib/intel64/ \
      -Wl,-rpath,/opt/intel/compilers_and_libraries_2018.0.128/linux/compiler/lib/intel64/
 
    Checks for "GenuineIntel" patched in the compiler and libraries.

    Minor failures in hydrogen_2P0_sfi.inp, hydrogen_1S_hhg_circular.inp, and
     hydrogen_1S_hhg_spline.inp

[6] AMD flang/clang version 11.0.0 (CLANG: AOCC_2.3.0-Build#85 2020_11_10), with the flags:
      -fopenmp -O3 -flto -march=native -funroll-loops -Mrecursive -mllvm -vector-library=LIBMVEC \
      -Kieee -fno-finite-math-only -cpp -D__BUILD_ID__='$(BUILD_ID)' \
      -Wl,-mllvm -Wl,-function-specialize -Wl,-mllvm -Wl,-region-vectorize \
      -Wl,-mllvm -Wl,-vector-library=LIBMVEC -Wl,-mllvm -Wl,-reduce-array-computations=3 \
      /opt/AMD/aocl/aocl-linux-aocc-2.2.0/lib/libflame.a \
      /opt/AMD/aocl/aocl-linux-aocc-2.2.0/lib/libblis.a \
      /opt/AMD/aocc-compiler-2.3.0/lib/libamdlibm.a \
      -lhugetlbfs -Wl,-z,max-page-size=2097152 -lmvec -ljemalloc -lm

   The compiler generates bad code for the timer routines; the implementation of system_clock()
   stops advancing real time after 2147.5 seconds; longer execution times are extracted from
   the time stamps reported by date_and_time intrinsic.

   Minor failures in hydrogen_1S_hhg_linear.inp, hydrogen_2P0_sfi.inp, hydrogen_1S_hhg_elliptical.inp,
   hydrogen_1S_hhg_circular.inp, hydrogen_1S_hhg_ell_refsol.inp, hydrogen_1S_hhg_spline.inp
