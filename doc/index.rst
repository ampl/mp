*mp* documentation
==================


Features
--------

* **Expression-based solver interface.**
  For solvers with an expression API,
  expression trees can be efficiently mapped. For example, AMPL expression
  ``max(a, b)`` directly maps to IBM ILOG Concert's ``IloMax``. The library
  has the following C++ interfaces of this kind, all of which support
  `AMPL extensions for logic and constraint programming`__:

  __ http://ampl.com/resources/logic-and-constraint-programming-extensions/

  - `Ilogcp <solvers/ilogcp>`_:
    IBM ILOG CPLEX and CPLEX CP Optimizer

  - `Gecode <solvers/gecode>`_

  - `JaCoP <solvers/jacop>`_

  - `LocalSolver <solvers/localsolver>`_

* **Conversion-based solver interface (WIP).**
  For solvers with more traditional 'flat' APIs, a customizable conversion
  layer translates expressions into constraints. For example, ``max(a, b)``
  is translated into ``<new var> = max(a, b)``, which is in turn redefined
  into a MIP construct or passed to the solver (Gurobi: ``GRBgenconstrMax``.)
  `Logical and CP constraints
  <http://ampl.com/resources/logic-and-constraint-programming-extensions/>`__
  are supported as well. For the solver API, an easy-to-adapt C++ wrapper class
  is provided. Currently it's two experimental interfaces:
  
  - `Gurobi <solvers/gurobidirect>`_ ('gurobidirect')

  - `IBM ILOG CPLEX <solvers/cplexdirect>`_ ('cplexdirect')

* **End-to-end solver testing script** for testing of various solver features.
  `Documentation. <test/end2end>`_

* An efficient type-safe C++ **adapter for the previous ASL library** for connecting solvers to AMPL and other systems:
  `source <src/asl>`_

* Reusable high-performance `.nl file reader <https://amplmp.readthedocs.io/en/latest/rst/nl-reader.html>`_
  which is up to `6x faster
  <http://zverovich.net/slides/2015-01-11-ics/socp-reformulation.html#/14>`_
  than the one provided by ASL. Documentation: https://amplmp.readthedocs.io/

* Database support on Linux and MacOS X.
  See `Database and spreadsheet connection guide`__.

  __  http://ampl.github.io/tables/

* `SMPSWriter <solvers/smpswriter>`_,
  a converter from deterministic equivalent of a two-stage stochastic
  programming (SP) problem written in AMPL to an SP problem in SMPS format.


Contents
========

.. toctree::
   :maxdepth: 2

   rst/nl-reader
   rst/problem
   rst/backend
   rst/reference



Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
