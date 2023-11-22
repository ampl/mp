Summary of recent updates to cplexmp for AMPL
============================================


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


