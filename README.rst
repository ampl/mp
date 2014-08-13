AMPL/MP
=======

.. image:: https://travis-ci.org/ampl/mp.png?branch=master
  :target: https://travis-ci.org/ampl/mp

.. image:: https://ci.appveyor.com/api/projects/status/91jw051om9q8pwt9
  :target: https://ci.appveyor.com/project/vitaut/mp

An open-source library for mathematical programming.

Features:

* Reusable high-performance `.nl <https://en.wikipedia.org/wiki/Nl_(format)>`_
  reader

* Efficient type-safe C++ API for connecting solvers to AMPL and other systems:
  `source <https://github.com/ampl/mp/tree/master/src/asl>`_

* Interfaces to solvers supporting
  `AMPL extensions for logic and constraint programming`__:

  __ http://ampl.com/resources/logic-and-constraint-programming-extensions/

  - IBM ILOG CPLEX and CPLEX CP Optimizer
    (`ilogcp <https://github.com/ampl/mp/tree/master/solvers/ilogcp>`_)

  - `Gecode <https://github.com/ampl/mp/tree/master/solvers/gecode>`_

  - `JaCoP <https://github.com/ampl/mp/tree/master/solvers/jacop>`_

* Interfaces to the following solvers:

  - `LocalSolver <https://github.com/ampl/mp/tree/master/solvers/localsolver>`_
  - `Sulum <https://github.com/ampl/mp/tree/master/solvers/sulum>`_

* Interfaces to other solvers via AMPL Solver Library

* Cross-platform build support with `CMake <http://www.cmake.org/>`_ and
  continuous integration systems. This includes third-party solvers and
  libraries (COIN-OR solvers with CMake support are available in the
  `ampl/coin <https://github.com/ampl/coin>`_ repository).

* `AMPLGSL <https://github.com/ampl/mp/tree/master/src/gsl>`_, an AMPL function
  library providing access to the GNU Scientific Library (GSL) functions.
  See the `AMPLGSL documentation <http://ampl.github.io/amplgsl>`_.

* Database support on Linux and MacOS X.
  See `Database and spreadsheet connection guide`__.

  __  http://ampl.github.io/tables/

* `SMPSWriter <https://github.com/ampl/mp/tree/master/solvers/smpswriter>`_,
  a converter from deterministic equivalent of a two-stage stochastic
  programming (SP) problem written in AMPL to an SP problem in SMPS format.

Binaries for the open-source AMPL solvers and libraries for major platforms
can be downloaded from the AMPL's `Open Source Solvers`__ page.

__ http://ampl.com/products/solvers/open-source/

Links
-----
`AMPL home <http://www.ampl.com/>`_ |
`AMPL book <http://ampl.github.io/ampl-book.pdf>`_ |
`Discussion group <https://groups.google.com/group/ampl>`_ |
`SolverStudio for Excel <http://solverstudio.org/languages/ampl/>`_

`AMPL models by Håkan Kjellerstrand <http://www.hakank.org/ampl/>`_
