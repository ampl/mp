AMPL/MP
=======

.. image:: https://travis-ci.org/ampl/mp.png?branch=master
  :target: https://travis-ci.org/ampl/mp

.. image:: https://ci.appveyor.com/api/projects/status/91jw051om9q8pwt9
  :target: https://ci.appveyor.com/project/vitaut/mp

An open-source library for mathematical programming.

Features:

* Reusable high-performance `.nl <https://en.wikipedia.org/wiki/Nl_(format)>`__ reader

* Efficient type-safe C++ API for connecting solvers to AMPL and other systems:
  `source <https://github.com/ampl/mp/tree/master/src/asl>`__

* Interfaces to solvers supporting
  `AMPL extensions for logic and constraint programming <http://ampl.com/resources/logic-and-constraint-programming-extensions/>`__:

  - IBM ILOG CPLEX and CPLEX CP Optimizer (`ilogcp <https://github.com/ampl/mp/tree/master/solvers/ilogcp>`__)

  - `Gecode <https://github.com/ampl/mp/tree/master/solvers/gecode>`__

  - `JaCoP <https://github.com/ampl/mp/tree/master/solvers/jacop>`__

* Interfaces to the following solvers:

  - `LocalSolver <https://github.com/ampl/mp/tree/master/solvers/localsolver>`__
  - `Sulum <https://github.com/ampl/mp/tree/master/solvers/sulum>`__

* Interfaces to other solvers via AMPL Solver Library

* Cross-platform build support with `CMake <http://www.cmake.org/>`__ and continuous integration
  systems. This includes third-party solvers and libraries (COIN-OR solvers with CMake support
  are available in the `ampl/coin <https://github.com/ampl/coin>`__ repository).

* `AMPLGSL <https://github.com/ampl/mp/tree/master/src/gsl>`__, an AMPL function
  library providing access to the GNU Scientific Library (GSL) functions.
  See the `AMPLGSL documentation <http://ampl.github.io/amplgsl>`__.

* Database support on Linux and MacOS X.
  See `Database and spreadsheet connection guide <http://ampl.github.io/tables/>`__.

* `SMPSWriter <https://github.com/ampl/mp/tree/master/solvers/smpswriter>`__, a converter
  from deterministic equivalent of a two-stage stochastic programming (SP) problem written in
  AMPL to an SP problem in SMPS format.

Binaries for the open-source AMPL solvers and libraries for major platforms can be downloaded
from the AMPL's `Open Source Solvers <http://ampl.com/products/solvers/open-source/>`__ page.

Links
-----
`AMPL home <http://www.ampl.com/>`__ | `AMPL book <http://ampl.github.io/ampl-book.pdf>`__ | `Discussion group <https://groups.google.com/group/ampl>`__ | `SolverStudio for Excel <http://solverstudio.org/languages/ampl/>`__ | `AMPL models by Håkan Kjellerstrand <http://www.hakank.org/ampl/>`__
