.. _solver-drivers:

Solver drivers
==============


.. _flat-solvers:

'Flat API' and new expression-based solvers
-----------------------------------------------

For solvers with traditional 'flat' (no expression trees) APIs,
non-linear AMPL expressions need to be reformulated.
For example, ``max(a, b)`` is translated into a constraint

.. code-block:: ampl

         new_var == max(a, b);

which is in turn reformulated for
MIP or passed to the solver natively (e.g., Gurobi: `GRBaddgenconstrMax`).

A late-2024 extension to MP allows passing expression trees
to solvers, see :ref:`supported-constraints`.

The implemented drivers are listed in :ref:`modeling-overview`.


.. _expression-solvers:

Old-API expression-based solvers
--------------------------------------

Some older drivers map directly from NL file's model, without
exploiting the automatic reformulation capabilites of MP.
For example, AMPL expression
``exp()`` maps to IBM ILOG Concert's ``IloExponent``.
The MP library
has the following C++ drivers of this kind, all of which support
`AMPL extensions for logic and constraint programming`__:

__ https://ampl.com/resources/logic-and-constraint-programming-extensions/

- `Ilogcp <https://github.com/ampl/mp/tree/develop/solvers/ilogcp>`_:
  IBM ILOG CPLEX and CPLEX CP Optimizer

- `Gecode <https://github.com/ampl/mp/tree/develop/solvers/gecode>`_

- `JaCoP <https://github.com/ampl/mp/tree/develop/solvers/jacop>`_ (unmaintained)

- `LocalSolver <https://github.com/ampl/mp/tree/develop/solvers/localsolver>`_ (unmaintained)

For constraint programming with MP's reformulations, use the MP-based
CP-SAT driver (see :ref:`modeling-overview`).


Specialized drivers
-------------------

- `SOCP solver <https://github.com/ampl/mp/tree/develop/solvers/cplexsocp>`_
  uses IBM ILOG CPLEX to solve problems convertable to SOCP form.

- `SSD solver <https://github.com/ampl/mp/tree/develop/solvers/ssdsolver>`_
  is a solver for problems with second-order stochastic dominance constraints.

- `SMPSWriter <https://github.com/ampl/mp/tree/develop/solvers/smpswriter>`_,
  a converter from deterministic equivalent of a two-stage stochastic
  programming (SP) problem written in AMPL to an SP problem in SMPS format.
