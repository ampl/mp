Summary of recent updates to MP2NL for AMPL
==========================================


## 20250814
- Changes in MP
  - Fixed reformulation of numberof and alldiff.
    Previously could produce an inefficient
    reformulation.
  - Option cvt:pre:boundlogarg (default 0) to bound
    arguments of logarithms nonnegative. Previously
    always done, sometimes deteriorating performance
    of nonlinear solvers.


## 20250128
- Updates for compilers compatibility


## 20241219
- Added *WantLogicalizedProd2Bin()*


## 20240926
- First release of the MP2NL driver
