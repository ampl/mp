Summary of recent updates to CP-SAT for AMPL
============================================

## unreleased
- Linear constraint coefficients must be integer: a fractional coefficient
  is an error, including a large one such as `1000000000.5`, which was
  rounded before. Values that differ from an integer only by a floating-point
  rounding error (within `1e-9`, e.g., `0.1*3*10`) are taken as that integer;
  this applies to coefficients, bounds, fixed values and exponents alike.
  Other fractional bounds of an integer-valued constraint are rounded inward.
  Finite constraint bounds of magnitude `1e18` and more are kept as given;
  those beyond the 64-bit integer range cannot be active and are accepted.
- An objective with a non-integer constant (e.g., `x + y + 2.5`) is
  accepted; before, the solve failed with "variable N is continuous".
- Fix: a solve after another solver's solve in the same AMPL session
  failed with "solution hint contains duplicate variables".
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
