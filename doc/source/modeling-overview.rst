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
   * - `Linear and mixed-integer <https://dev.ampl.com/solvers/index.html#linear-and-mixed-integer-solvers>`__
     - `Gurobi <https://dev.ampl.com/solvers/gurobi/index.html#features>`__, `CPLEX <https://dev.ampl.com/solvers/cplex/index.html#features>`__, `Xpress <https://dev.ampl.com/solvers/xpress/index.html#features>`__, `COPT <https://dev.ampl.com/solvers/copt/index.html#features>`__,
       `HiGHS <https://dev.ampl.com/solvers/highs/index.html#features>`__, `SCIP <https://dev.ampl.com/solvers/scip/index.html#features>`__, `MOSEK <https://dev.ampl.com/solvers/mosek/index.html#features>`__, `CBC <https://dev.ampl.com/solvers/cbc/index.html#features>`__,
       `GCG <https://dev.ampl.com/solvers/gcg/index.html#features>`__, `cuOpt <https://dev.ampl.com/solvers/cuopt/index.html#features>`__
     - Logic becomes indicator or general constraints where the solver has
       them, big-M rows otherwise; other nonlinear functions are
       approximated piecewise-linearly.
   * - Conic
     - MOSEK, COPT, Gurobi, Xpress, CPLEX, SCIP
     - Recognizes second-order and exponential cones from their algebraic
       form (:ref:`socp`).
   * - `Global nonlinear <https://dev.ampl.com/solvers/index.html#nonlinear-solvers>`__
     - Gurobi, Xpress, SCIP, `BARON <https://dev.ampl.com/solvers/baron/index.html#features>`__ (``baronmp``);
       COPT solves nonlinear models locally (``acc:_expr=1``)
     - Passes nonlinear functions as expression trees
       (:ref:`flat-vs-expressions`).
   * - `Constraint programming <https://dev.ampl.com/solvers/index.html#constraint-programming-solvers>`__
     - `CP-SAT <https://dev.ampl.com/solvers/cpsat/index.html#features>`__ (integer models)
     - Receives logic, counting, ``alldiff``, products, ``div`` and ``mod``
       natively. The same AMPL constructs are reformulated for MIP solvers,
       so a constraint-programming model can also be solved by any MIP solver.
   * - Any NL solver
     - :ref:`MP2NL <mp2nl>` (e.g. Knitro, Ipopt, Couenne, BARON)
     - Applies MP's reformulations, then calls the NL solver.

Each solver's page on `AMPL Development <https://dev.ampl.com/solvers/index.html>`_
(linked above) summarizes the problem types and expressions it receives natively,
its algorithms and options. The exact list for an installed driver is in its
:ref:`acceptance options <acceptance-options>` (``<solver> -=acc:``); other
features are compared in :ref:`support-by-solvers`. Binaries are available
through the `AMPL Portal <https://portal.ampl.com>`_.

**Where to read next:**
:ref:`modeling-formulations` (choosing between formulations),
:ref:`efficiency` and :ref:`numerical_accuracy` (writing models that solve
well), :ref:`modeling-tools` (reformulation options, exploring what the
solver receives, solution check, MP2NL), and :ref:`modeling-troublesh`.
Small modeling tasks are covered in the
`MP Modeling Series <https://ampl.com/streamlit/Modeling_Tips>`_.
