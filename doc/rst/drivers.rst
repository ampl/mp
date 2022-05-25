.. _solver_drivers:

Solver drivers
==============

Expression-based APIs
---------------------

Expression-based solver APIs can efficiently map
NL forests.
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


Flat APIs
---------

For solvers with traditional "flat" APIs, non-linear AMPL expressions need
to be reformulated.
For example, ``max(a, b)`` is translated into a constraint meaning
``<new var> = max(a, b)``, which is in turn reformulated for
MIP or passed to the solver natively (Gurobi: `GRBaddgenconstrMax`).
See the :ref:`modeling_guide`.

Currently there are three experimental implementations:

- `x-gurobi <https://github.com/ampl/mp/tree/master/solvers/gurobidirect>`_
  (available in the AMPL distribution bundle)

- `x-cplex <https://github.com/ampl/mp/tree/master/solvers/cplexdirect>`_

- `copt <https://github.com/ampl/mp/tree/master/solvers/copt>`_
  (available in the AMPL distribution bundle)

- `HiGHS <https://github.com/ampl/mp/tree/develop/solvers/highsdirect>`_
  (see the `HiGHS website <https://highs.dev/>`_)


Specialized drivers
-------------------

- `SOCP solver <https://github.com/ampl/mp/tree/master/solvers/cplex>`_
  uses IBM ILOG CPLEX to solve problems convertable to SOCP form.

- `SSD solver <https://github.com/ampl/mp/tree/master/solvers/ssdsolver>`_
  is a solver for problems with second-order stochastic dominance constraints.

- `SMPSWriter <https://github.com/ampl/mp/tree/master/solvers/smpswriter>`_,
  a converter from deterministic equivalent of a two-stage stochastic
  programming (SP) problem written in AMPL to an SP problem in SMPS format.
