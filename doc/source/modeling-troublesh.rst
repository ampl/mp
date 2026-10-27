
.. _modeling-troublesh:

Troubleshooting
---------------------------

Solving efficiency
*********************

Always look for a more efficient model.
Bounds could be tightened; coefficients' magnitude can be reduced
by rescaling of the data; you might not need certain constraints
or logical conditions when using a different approach.
Sometimes, just a different
formulation can help, for example, a logical condition
can be manually linearized.
See sections :ref:`efficiency` and :ref:`numerical_accuracy`; for choosing
between whole formulations, and checking a candidate before a benchmark, see
:ref:`modeling-formulations` and :ref:`before-benchmarking`.

Play with `AMPL options <https://dev.ampl.com/ampl/options.html>`_
and :ref:`solver-options`.
Prominent ones are AMPL presolve
(switch off: ``ampl: option presolve 0;``), solver's presolve
(``ampl: option gurobi_options 'presolve=0';``) and others
(tolerances, *numfocus*, *intfocus*, automatic tuning, etc.), as well as
MP reformulation options, see :ref:`supported-constraints`.


Exporting the solver model
*********************************

To see what MP and/or the solver do with your model, export
the solver's received model, and, if possible, the solver's presolved model:

.. tab-set::

   .. tab-item:: Python

      .. code-block:: python

         solver = "gurobi"                          # any MP solver
         ampl.option[solver + "_auxfiles"] = "rc"   # pass AMPL names to the solver
         ampl.solve(solver=solver, mp_options="writeprob=disj.lp")

   .. tab-item:: AMPL

      .. code-block:: ampl

         ampl: option solver gurobi;           ## any MP solver
         ampl: option gurobi_auxfiles rc;      ## pass AMPL names to the solver
         ampl: option mp_options 'writeprob=disj.lp';
         ampl: solve;

Names appear in the file where the solver's export format keeps them (the
HiGHS LP export, for example, uses generic names). Solvers that can also export
their presolved model accept ``writepresolved`` (for example
``gurobi_options="writeprob=disj.lp writepresolved=disj_pre.lp"``).

MP offers ways to :ref:`explore automatic reformulations<explore-reformulations>`.
Formats and options are in :ref:`writeprob`; how to read the big-M constants
MP generated is shown in :ref:`declared-bounds-bigm`.


.. _reading-solver-results:

Reading solver results
*************************

A solver status is a claim. The status values themselves (``solve_result``,
``solve_result_num`` and their ranges) are listed in :ref:`solve-result-codes`.
Check a status before acting on it:

- **"Optimal" from a local solver is a local claim.** On a nonconvex model, a
  local NLP or MINLP solver's ``solved`` means a point satisfying local
  optimality conditions, not a proven optimum. On linear, convex quadratic,
  conic and MILP models the same status is a proof.
- **A time-limit stop can return the starting point.** If a solver stops
  without an improved point, the variables keep their previous values and
  look like a solution. Compare them with the values before the solve.
- **"Infeasible" can be wrong** on hard nonconvex or numerically difficult
  models. Before acting on it, solve with a second solver, or solve a version
  with the constraints made elastic (slack variables with a large penalty)
  and look at which slacks are needed. MP solvers can build that version
  themselves (:ref:`feasibilityrelaxation`) or report an irreducible
  infeasible subset (:ref:`iis`).
- **Check the point itself**: the :ref:`automatic solution check <solution-check>`
  warnings, and the violated logical constraints
  (``display {i in 1.._nlogcons: !_logcon[i].val} (_logconname[i]);``).
- **Do not screen results with the reported gap alone.** Recompute the gap
  from the objective and the bound where it matters; see
  :ref:`returnMIPgap` and :ref:`returnBestBound` for retrieving them.
- **Cross-check bounds against incumbents.** When several runs (solvers,
  seeds, formulations) solve the same instance, every feasible objective is an
  incumbent and every dual bound bounds the same optimum, so an incumbent
  better than any run's bound exposes a defect in one of the two runs.
  Ignore differences below about :math:`10^{-9}` relative, and those within
  the solvers' optimality tolerance.

In a script, an unsuccessful solve may not raise an exception: check
``solve_result`` after every solve.

.. tab-set::

   .. tab-item:: Python

      .. code-block:: python

         before = ampl.var["x"].to_dict()
         ampl.solve(solver="highs", highs_options="lim:time=60")
         after = ampl.var["x"].to_dict()
         moved = any(abs(after[k] - before[k]) > 1e-9 for k in after)
         if ampl.solve_result not in ("solved", "limit") or not moved:
             print(ampl.get_value("solve_message"))

   .. tab-item:: AMPL

      .. code-block:: ampl

         solve;
         if solve_result not in {'solved', 'limit'} then
            print solve_message;


.. _comparing-formulations:

Comparing formulations and solvers
**************************************

How to choose and check the formulations themselves (root bounds, indicators
or big-M, known reformulations) is in :ref:`modeling-formulations`; this
section is about measuring them fairly.

- **Validate first.** A reformulation is trusted only after it reproduces
  the original's objective and status on identical data, at more than one
  instance size and with more than one solver.
- **Repeat with several random seeds** (option *tech:seed* where available).
  Mixed-integer solving times vary a lot under changes that should not matter
  (Lodi & Tramontani 2013); a 2.4 times speedup seen on one seed disappeared
  over five.
- **Compare on the same instances.** Compare times only over instances both
  sides solved, and report how many there were. A mean over each side's own
  solved instances rewards a formulation for failing on the hard ones; in one
  benchmark it turned "1.4 times slower" into "1.6 times faster".
- **Keep the axes apart**: number of instances solved, final gap, and time to
  optimality often disagree.
- **Summarize robustly**: shifted geometric means of times, or performance
  profiles (Dolan & Moré 2002).
- **Run every arm under the same conditions**: thread count, time limit,
  options, starting values (:ref:`warm-start`), and machine load. Record the solver and driver versions with every
  result.
- **Change one thing per arm.** If an arm changes two things that may both
  matter, add the missing combination rather than argue which one did it.

*References:* E.D. Dolan, J.J. Moré (2002). Benchmarking optimization software
with performance profiles. *Math. Programming* 91, 201-213. A. Lodi,
A. Tramontani (2013). Performance variability in mixed-integer programming.
*INFORMS TutORials in Operations Research*, doi:10.1287/educ.2013.0112.


Contacting support
*******************************

If you decide to contact AMPL or solver support, please provide a (possibly reduced)
version of your model reproducing the issue. Please also provide
the AMPL solver version by running ``highs -v``
or solver logs obtained with options ``version outlev=1 timing=2``.
From a Python script, write the model and data of the failing case from a
**fresh** AMPL session (a session in which the solver crashed may not respond):

.. code-block:: python

    ampl = AMPL()
    ampl.read("model.mod")
    # ... load the same data as in the failing run ...
    ampl.export_model("repro.mod")
    ampl.export_data("repro.dat")

If you cannot reduce the model and don't
want to show it, give us just the NL file produced by

.. tab-set::

   .. tab-item:: Python

      .. code-block:: python

         ampl.write("bdisj")      # binary NL file disj.nl, without auxiliary files

   .. tab-item:: AMPL

      .. code-block:: ampl

         ampl: option auxfiles ''; write bdisj;

For solver support, use MPS or another full-precision model format.
