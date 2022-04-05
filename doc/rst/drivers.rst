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
Currently there are two experimental implementations:

- `x-gurobi <https://github.com/ampl/mp/tree/master/solvers/gurobidirect>`_
  (available in the AMPL distribution bundle; see the :ref:`modeling_guide`)

- `x-cplex <https://github.com/ampl/mp/tree/master/solvers/cplexdirect>`_

