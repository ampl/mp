Reusable components and solver driver API
=========================================


NL file reader
--------------

High-performance :doc:`.nl file reader <nl-reader>`
which is up to `6x faster
<http://zverovich.net/slides/2015-01-11-ics/socp-reformulation.html#/14>`_
than the one provided by the traditional
`AMPL Solver Library (ASL)
<https://ampl.com/resources/learn-more/hooking-your-solver-to-ampl/>`_.
It can be customized for most efficient translation of NL format into
solver API using the `mp::ProblemBuilder` concept.

* Alternatively, convenience classes `mp::Problem` and `mp::ColProblem`
  can be used for intermediate storage of the NL model. From there,
  `mp::ExprVisitor` and `mp::ModelFlattener` walk the NL forest.

:doc:`More info on NL Reader <nl-reader>`


NL model management
-------------------

* Class `mp::BasicModelManager` standardizes model input and results output.


Model conversion and presolve
-----------------------------

* Classes `mp::FlatConverter` and `mp::MIPFlatConverter` facilitate conversion of
  NL expressions which are not natively accepted by a solver into simpler forms.

  * `Logical and CP constraints
    <http://ampl.com/resources/logic-and-constraint-programming-extensions/>`__
    are supported.

* Class `mp::pre::Presolver` pre- and postsolves solutions and suffixes.

:doc:`More info on model / solution conversions <flatcvt>`


Solver logic
------------

* Classes `mp::BasicSolver`, `mp::Backend` and `mp::MIPBackend`
  standardize solver behavior such as common options and suffixes
  and are recommended for new interfaces.

:doc:`More info on solver logic API <backend>`


API reference
-------------

:doc:`Reference <reference>`

#.. include:: nl-reader.rst
#.. include:: flatcvt.rst
#.. include:: backend.rst
#.. include:: reference.rst
