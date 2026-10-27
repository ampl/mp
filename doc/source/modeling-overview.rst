.. _modeling-overview:

Modeling overview and implemented solvers
------------------------------------------

Solvers built with AMPL's MP library (`ampl/mp <https://github.com/ampl/mp>`_)
accept a much wider range of model expressions than a solver's own API:
you write the model the way you think about it, and MP passes each construct
to the solver natively where the solver supports it, and reformulates it
otherwise.

**What MP handles**, with details in :ref:`expressions_supported`:

- conditional and logical constraints: ``if-then-else``, ``==>``, ``or``,
  ``and``, ``exists``, ``forall``;
- counting and constraint-programming constructs: ``count``, ``atmost``,
  ``atleast``, ``exactly``, ``numberof``, ``alldiff``, ``!=``;
- piecewise-linear functions, ``abs``, ``min``, ``max``;
- rounding and integer division: ``floor``, ``ceil``, ``round``, ``trunc``,
  ``div``, ``mod``;
- nonlinear functions and products, second-order and exponential cones;
- complementarity (``complements``) and set membership (``in``; being
  reworked, see :ref:`set-membership`);
- multiple objectives (:ref:`multiple-objectives`).

**Solver classes.** What happens to a construct depends on the solver:

.. list-table::
   :header-rows: 1
   :widths: 22 38 40

   * - Class
     - MP-based drivers
     - What MP does
   * - Linear and mixed-integer
     - Gurobi, CPLEX, Xpress, COPT, HiGHS, SCIP, MOSEK, CBC, GCG, cuOpt
     - Logic becomes indicator or general constraints where the solver has
       them, big-M rows otherwise; other nonlinear functions are
       approximated piecewise-linearly.
   * - Conic
     - MOSEK, COPT, Gurobi, Xpress, CPLEX, SCIP
     - Recognizes second-order and exponential cones from their algebraic
       form (:ref:`socp`).
   * - Global nonlinear
     - Gurobi, Xpress, SCIP, BARON (``baronmp``); COPT solves nonlinear
       models locally (``acc:_expr=1``)
     - Passes nonlinear functions as expression trees
       (:ref:`flat-vs-expressions`).
   * - Constraint programming
     - CP-SAT (integer models)
     - Receives logic, counting, ``alldiff``, products, ``div`` and ``mod``
       natively. The same AMPL constructs are reformulated for MIP solvers,
       so a constraint-programming model can also be solved by any MIP solver.
   * - Any NL solver
     - :ref:`MP2NL <mp2nl>` (e.g. Knitro, Ipopt, Couenne, BARON)
     - Applies MP's reformulations, then calls the NL solver.

Which constructs a particular solver accepts natively is listed in its
:ref:`acceptance options <acceptance-options>` (``<solver> -=``); its other
features are summarized in :ref:`support-by-solvers`. Binaries are available
through the `AMPL Portal <https://portal.ampl.com>`_, and each solver's page
on `AMPL Development <https://dev.ampl.com/solvers/index.html>`_ describes its
options and solver-specific features.

**Where to read next:**
:ref:`modeling-formulations` (choosing between formulations),
:ref:`efficiency` and :ref:`numerical_accuracy` (writing models that solve
well), :ref:`modeling-tools` (reformulation options, exploring what the
solver receives, solution check, MP2NL), and :ref:`modeling-troublesh`.
Small modeling tasks are covered in the
`MP Modeling Series <https://ampl.com/streamlit/Modeling_Tips>`_.
