Summary of recent updates to CP-SAT for AMPL
============================================

## unreleased
- Option *tech:varbound* is now called *cvt:mip_max_bound*;
  *tech:varbound* and *varbound* remain accepted.
- Linear constraint coefficients must be integer: a fractional coefficient
  is an error, including a large one such as `1000000000.5`, which was
  rounded before. A value that differs from an integer only by
  floating-point rounding, such as `0.1*3*10`, is taken as that integer:
  within `1e-9`, or for magnitudes above about `1e6` within 4 units in the
  last place of the value, but never by `0.25` or more (a value ending in
  `.5` is never taken as an integer). This applies to coefficients, bounds
  and fixed values alike. Other fractional bounds of an integer-valued
  constraint are rounded inward. Before, a coefficient within a relative
  `1e-9` of an integer was rounded (`1e7 + 0.001` became `1e7`); such a
  coefficient is now an error. Use `round()` in the model if rounding is
  intended.
- Finite constraint bounds of magnitude `1e18` and more are kept as given.
  A bound beyond the 64-bit integer range is no longer an error: an upper
  bound above the range (a lower bound below it) never binds; a lower bound
  above the range (an upper bound below it) makes the model infeasible.
- Fix: *cvt:mip_max_bound* (default `1e12`) now replaces infinite variable
  bounds only; finite bounds are passed as given. Before, a finite bound
  of that magnitude or more was replaced too: an objective constant such
  as in `x + y + 2e12`, or a lower bound such as `x >= 2e12`, made the
  model infeasible, and an upper bound such as `x <= 2e14` was silently
  lowered to `1e12`. An infinite bound opposite a finite one beyond
  *cvt:mip_max_bound* is set *cvt:mip_max_bound* beyond it.
- Warnings: "VARBOUND" lists the infinite bounds replaced by
  *cvt:mip_max_bound* and recommends finite bounds for all variables;
  "VARBOUND_ACTIVE" lists solution values at such a bound, where the
  result may depend on *cvt:mip_max_bound* (e.g., an unbounded model).
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
