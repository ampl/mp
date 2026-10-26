Summary of recent updates to CP-SAT for AMPL
============================================

## 20261101
- First release of the CP-SAT driver.
- Products and integer powers of integer variables are passed
  natively (*acc:mul*, *acc:powc*); before, they were approximated
  with continuous variables, which CP-SAT rejects.
- *div* and *mod* are passed natively (*acc:intdiv*, *acc:mod*).
- Initial values of the variables (including MP's auxiliary
  variables) are passed as CP-SAT solution hints (*alg:start*).
