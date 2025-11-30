Summary of recent updates to MP2NL for AMPL
==========================================


## 20251201
- Change in MP: option *obj:multi:options*: fix integer-valued
  objective-specific options.


## 20251121
- New option *obj:multi:options* to control
  whether multi-objective option suffixes
  are used.


## 20251021
- Changes in MP: option *cvt:expr:nlassign*.


## 20250814
- Changes in MP
  - Improved preprocessing of logical
    and combinatorial expressions
    (options cvt:pre:unnest, cvt:pre:sort).
  - Option cvt:pre:boundlogarg (default 0) to bound
    arguments of logarithm nonnegative. Previously
    always done, sometimes deteriorating performance
    of nonlinear solvers.


## 20250128
- Updates for compilers compatibility


## 20241219
- Added *WantLogicalizedProd2Bin()*


## 20240926
- First release of the MP2NL driver
