AMPL/MP
=======

An open-source library for mathematical programming.

Components:

* Reusable high-performance `.nl <https://en.wikipedia.org/wiki/Nl_(format)>`__ reader

* Efficient type-safe C++ API for writing AMPL solvers:
  `source <https://github.com/ampl/mp/tree/master/src/asl>`__

* AMPL bindings for the GNU Scientific Library: `docs <http://ampl.github.io/amplgsl/>`__,
  `source <https://github.com/ampl/mp/tree/master/src/gsl>`__

* Interfaces to solvers supporting
  `AMPL extensions for logic and constraint programming <http://ampl.com/resources/logic-and-constraint-programming-extensions/>`__:

  - IBM ILOG CPLEX and CPLEX CP Optimizer (`ilogcp <https://github.com/ampl/mp/tree/master/solvers/ilogcp>`__)

  - `Gecode <https://github.com/ampl/mp/tree/master/solvers/gecode>`__

  - `JaCoP <https://github.com/ampl/mp/tree/master/solvers/jacop>`__

* Interfaces to the following solvers:

  - `LocalSolver <https://github.com/ampl/mp/tree/master/solvers/localsolver>`__
  - `Sulum <https://github.com/ampl/mp/tree/master/solvers/sulum>`__

* Interfaces to other solvers via AMPL Solver Library

Links
-----
`AMPL home <http://www.ampl.com/>`__ | `AMPL book <http://ampl.github.io/ampl-book.pdf>`__ | `Discussion group <https://groups.google.com/group/ampl>`__ | `SolverStudio for Excel <http://solverstudio.org/languages/ampl/>`__ | `AMPL models by Håkan Kjellerstrand <http://www.hakank.org/ampl/>`__