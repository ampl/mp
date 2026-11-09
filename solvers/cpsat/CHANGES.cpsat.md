Summary of recent updates to CP-SAT for AMPL
============================================

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
