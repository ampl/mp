*mp* documentation
==================

MP Library is a set of tools recommended to create new AMPL solver interfaces.

Features
--------

Reusable building blocks for new interfaces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* High-performance :doc:`.nl file reader <rst/nl-reader>`
  which is up to `6x faster
  <http://zverovich.net/slides/2015-01-11-ics/socp-reformulation.html#/14>`_
  than the one provided by ASL. It can be customized for most efficient translation of NL format into
  solver API.

* Classes `mp::Backend` and `mp::MIPBackend`
  standardize solver behavior such as common options and suffixes
  and are recommended for new interfaces.

* Classes `mp::FlatConverter` and `mp::MIPFlatConverter` facilitate conversion of
  NL expressions which are not natively accepted by a solver into simpler forms.
  `Logical and CP constraints
  <http://ampl.com/resources/logic-and-constraint-programming-extensions/>`__
  are supported.

* Convenience classes `mp::Problem` and `mp::ColProblem` can be used for
  intermediate storage of the NL model.
  `mp::ExprVisitor` and `mp::ExprFlattener` walk NL forest top-down.

Concrete solver interfaces
^^^^^^^^^^^^^^^^^^^^^^^^^^

* Interfaces to solvers with **expression-based APIs.**
  For solvers with an expression-based API,
  NL forests can be efficiently mapped. For example, AMPL expression
  ``exp()`` maps to IBM ILOG Concert's ``IloExponent``. The library
  has the following C++ interfaces of this kind, all of which support
  `AMPL extensions for logic and constraint programming`__:

  __ http://ampl.com/resources/logic-and-constraint-programming-extensions/

  - `Ilogcp <https://github.com/ampl/mp/tree/master/solvers/ilogcp>`_:
    IBM ILOG CPLEX and CPLEX CP Optimizer

  - `Gecode <https://github.com/ampl/mp/tree/master/solvers/gecode>`_

  - `JaCoP <https://github.com/ampl/mp/tree/master/solvers/jacop>`_

  - `LocalSolver <https://github.com/ampl/mp/tree/master/solvers/localsolver>`_

* Interfaces to solvers with **"flat" APIs** (WIP).
  For solvers with more traditional "flat" APIs, class `mp::MIPFlatConverter`
  reformulates many non-linear AMPL expressions.
  For example, ``max(a, b)`` is translated into a constraint meaning
  ``<new var> = max(a, b)``, which is in turn reformulated for
  MIP or passed to the solver natively (Gurobi: `GRBaddgenconstrMax`).
  Currently there are two experimental implementations:
  
  - `Gurobi <https://github.com/ampl/mp/tree/master/solvers/gurobidirect>`_

  - `IBM ILOG CPLEX <https://github.com/ampl/mp/tree/master/solvers/cplexdirect>`_

Other utilities
^^^^^^^^^^^^^^^

* An efficient type-safe C++ **adapter for the previous ASL library** for
  connecting solvers to AMPL and other systems:
  `source <../../../src/asl>`_

* `SMPSWriter <../../../solvers/smpswriter>`_,
  a converter from deterministic equivalent of a two-stage stochastic
  programming (SP) problem written in AMPL to an SP problem in SMPS format.

* **End-to-end solver testing script** for testing of various solver features:
  `source <../../../test/end2end>`_


Contents
--------

.. toctree::
   :maxdepth: 2

   rst/nl-reader
   rst/backend
   rst/flatcvt
   rst/reference



Indices and tables
------------------

* :ref:`genindex`
* :ref:`search`
