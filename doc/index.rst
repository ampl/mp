*mp* documentation
==================

MP Library is a set of solver drivers and tools
recommended to create new AMPL solver drivers.

Features
--------

Solver drivers and a modeling guide
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Drivers for solvers with **expression-based APIs.**
  In this case,
  NL forests can be efficiently mapped to the solver API.
  For example, AMPL expression
  ``exp()`` maps to IBM ILOG Concert's ``IloExponent``. The library
  has the following C++ drivers of this kind, all of which support
  `AMPL extensions for logic and constraint programming`__:

  __ http://ampl.com/resources/logic-and-constraint-programming-extensions/

  - `Ilogcp <https://github.com/ampl/mp/tree/master/solvers/ilogcp>`_:
    IBM ILOG CPLEX and CPLEX CP Optimizer

  - `Gecode <https://github.com/ampl/mp/tree/master/solvers/gecode>`_

  - `JaCoP <https://github.com/ampl/mp/tree/master/solvers/jacop>`_

  - `LocalSolver <https://github.com/ampl/mp/tree/master/solvers/localsolver>`_

* Drivers for solvers with traditional **"flat" APIs**.
  For solvers with "flat" APIs, non-linear AMPL expressions need
  to be reformulated.
  For example, ``max(a, b)`` is translated into a constraint meaning
  ``<new var> = max(a, b)``, which is in turn reformulated for
  MIP or passed to the solver natively (Gurobi: `GRBaddgenconstrMax`).
  Currently there are two experimental implementations:

  - `x-gurobi <https://github.com/ampl/mp/tree/master/solvers/gurobidirect>`_
    (available in the AMPL distribution bundle)

  - `x-cplex <https://github.com/ampl/mp/tree/master/solvers/cplexdirect>`_

* A :doc:`summary guide for efficient modeling for x-gurobi and x-cplex <rst/model-guide>`.


Reusable building blocks for new interfaces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* High-performance :doc:`.nl file reader <rst/nl-reader>`
  which is up to `6x faster
  <http://zverovich.net/slides/2015-01-11-ics/socp-reformulation.html#/14>`_
  than the one provided by ASL. It can be customized for most efficient translation of NL format into
  solver API.

* Classes `mp::FlatConverter` and `mp::MIPFlatConverter` facilitate conversion of
  NL expressions which are not natively accepted by a solver into simpler forms.

  * `Logical and CP constraints
    <http://ampl.com/resources/logic-and-constraint-programming-extensions/>`__
    are supported.

* Classes `mp::BasicSolver`, `mp::Backend` and `mp::MIPBackend`
  standardize solver behavior such as common options and suffixes
  and are recommended for new interfaces.

* Convenience classes `mp::Problem` and `mp::ColProblem` can be used for
  intermediate storage of the NL model.
  From there, `mp::ExprVisitor` and `mp::ExprFlattener` walk NL forest top-down.

Other utilities
^^^^^^^^^^^^^^^

* An efficient type-safe `C++ adapter for the previous ASL library
  <https://github.com/ampl/mp/tree/master/src/asl>`_ for
  connecting solvers to AMPL and other systems.

* `SMPSWriter <https://github.com/ampl/mp/tree/master/solvers/smpswriter>`_,
  a converter from deterministic equivalent of a two-stage stochastic
  programming (SP) problem written in AMPL to an SP problem in SMPS format.

* `End-to-end solver testing script
  <https://github.com/ampl/mp/tree/master/test/end2end>`_ for testing of
  various solver features.


Contents
--------

.. toctree::
   :maxdepth: 2

   rst/model-guide
   rst/nl-reader
   rst/flatcvt
   rst/backend
   rst/reference
   rst/troublesh


Indices and tables
------------------

* :ref:`genindex`
* :ref:`search`
