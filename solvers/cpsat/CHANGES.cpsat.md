Summary of recent updates to CP-SAT for AMPL
============================================

## unreleased
- Fixed integer variables are kept integer (*cvt:pre:continuous_fixed_vars*
  is 0 by default for CP-SAT).
- Models with contradictory variable bounds, or with an integer variable
  fixed at a fractional value, are reported infeasible instead of failing
  with an error.

## 20261029
- First release of the CP-SAT driver.
- Exact products and integer powers of integer variables are supported
  natively (*acc:mul*, *acc:powc*).
- *div* and *mod* are passed natively (*acc:intdiv*, *acc:mod*).
- Initial variable values are passed as CP-SAT solution hints
  (*alg:start*).
- New option *tech:logfile*: appends the CP-SAT search log to a file,
  regardless of *tech:outlev* (which only controls whether the log is
  also written to stdout). Without it, behavior is unchanged.
