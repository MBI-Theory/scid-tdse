Last updated: 2022 Oct 08
------------

SCID-TDSE mostly works on macOS, achieving excellent performance on Apple
M-series CPUs. There are however a number of factors to consider, namely:

* You will need a working Fortran compiler and build tools. Homebrew
(https://brew.sh/) is a convenient way to install those; just follow
instructions on the web site. Once homebrew is installed, you will need
to install a few formulae/casks:

  brew install gcc openblas xquartz

(openblas and xquartz are optional, but recommended)

If you wish to use OpenDX visualization, you can access it through an Unbuntu
guest VM, which can be installed with "brew install multipass".

* macOS has relatively low maximum (and even lower default) stack size, which
may be insufficient to run large calculations. All examples in the test set
run successfully with a 3-megabyte stack (main thread and OpenMP), but larger
stacks may be needed for other inputs.

* The GNU compilers are still not fully stable on Apple M1/M2 CPUs. gfortran
12.2 _mostly_ works, but may still experience catastrophic test failures at
high optimization levels. In particular, all optimization levels beyond -O1
produce bad code, which experiences segmentation faults for large inputs.  An
example configuration file can be found in "macos_m1-gfortran_opt.mak" in the
configs/ subdirectory.  This configuration should also work on Intel Macs - but
this was not tested.

* SCID-TDSE compiled with gfortran on M1 macs is still not completely stable;
crashes do happen for some input and some test cases, including
hydrogen_1S_hhg_spline, hydrogen_1S_hhg_elliptical, mhydrogen_1S_hhg_circular,
hydrogen_1S_hhg_ell_refsol. The problems appear to be probabilistic in nature,
suggesting a race condition. If your calculation is affected by the crashes,
try reducing the number of OpenMP threads, or disabling OpenMP altogether.

* SCID-TDSE requires BLAS and LAPACK. The versions provided by the Accelerate
framework work, and and are configured by default in macos_m1-gfortran_opt.mak.
openblas also works, and provides somewhat better performance.

* On M1/M2 Macs, not all cores are equal. Requesting a number of cores larger
than the number of available performance cores may lead to erratic performance
results.

* Due to the aggressive macOS power management, it is recommended to start
calculations using "caffeinate" wrapper. Otherwise, inconsistent timing and
poor parallel performance are likely.

