Summary of recent updates to CPLEX MP for AMPL
==============================================

## TBD
- Added fixed model for `mip:basis`


## 20231122
- Public beta; implemented most [features](https://mp.ampl.com/features-guide.html) in MP:
  model export, warm start, input and output basis, feasibility relaxation,
  multiple solutions, kappa, unbounded rays, IIS, return MIP gap, 
  return dual bound
- Native Model support:
    - continuous and integer variables
    - multiple objectives, quadratic objectives
    - linear and quadratic constraints, indicator constraints and PL functions
    - special ordered sets (type 1 and 2)


## 20231117
- MP update: fixed graceful exit on Ctrl-C from AMPL in Linux
  and fixed issue with reading text-format NL files


## 20230815
- Fixed a bug causing repeated names for
  auxiliary variables and constraints.
- Option values can be assigned without '='.
- Changed default tolerance for strict comparisons
  to 0 (option cmp:eps, #102.)
- Fixed a bug where equivalent conditional
  comparisons were not unified.

 