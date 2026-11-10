.. _multiple-objectives:


Multiple objectives
----------------------------------

.. image:: images/berg-tal.svg
  :width: 200
  :align: right

Many real-world problems have multiple objectives; often this scenario
is tackled by blending all the objectives
by linear combination when formulating the model, or by minimizing
unwanted objective deviations from a pre-specified goal.

Alternatively, if some constraints are 'soft', they can be modeled
by penalties included in the primary or secondary objectives.

To consider multiple objectives in an AMPL model,
either natively supported by the solver, or emulated, use
solver option :ref:`obj:multi (multiobj) <multipleObjectives>`.
Otherwise, only the 1st objective is considered
(or any objective specified by ``obj:no``.)

.. code-block:: ampl

    minimize total_cost {s in STORE}:
       sum {j in FOOD} cost[s,j] * Buy[j];

    minimize total_number:  sum {j in FOOD} Buy[j];

Objectives can be blended, or hierarchical (lexicographical).


Blended objectives
************************************************

By default, all objectives are blended together
(summed up with the corresponding signs.)
`Suffixes <https://dev.ampl.com/ampl/reference/suffixes.html>`__
``.objweight`` can be used to change the individual weights
and objective senses, according to the option ``obj:multi:weight``.
Example:

.. code-block:: ampl

    suffix objweight IN;                       # Declare the suffix
    let {s in STORE} total_cost.objweight := TCFactor[s];
    let total_number.objweight := 1;


Lexicographical objectives
********************************************************

To apply hierarchical optimization, use suffix ``.objpriority``,
as described in the ``obj:multi`` option description.

First, declare the suffixes:

.. code-block:: ampl

    suffix objpriority IN, integer;
    suffix objabstol IN;
    suffix objreltol IN;

    maximize ReverseSeniority {e in 1..2, i in I: E[i]==e}:
      sum {t in V[i]: Pr[i, t]==0}
        S[i] * x[i, t]
      suffix objpriority (2-e)*S_range + 1 + S[i] - min {j in I} S[j];

Suffixes ``.objabstol`` and ``.objreltol`` allow for objective degradation.
However their exact meaning can vary for a solver's native multi-objective
mode (``obj:multi=1``), in particular for LPs. Consult
:ref:`Feature support by solvers <support-by-solvers>`
and
`AMPL Colab notebooks with lexicographical objectives <https://colab.ampl.com/tags/lexicographic-objectives.html>`__.

- **A tolerance of 0 means "not set".** AMPL suffixes have default value 0,
  so ``.objreltol 0`` leaves the solver's default degradation in place,
  which need not be 0. To allow no degradation, set a very small value,
  such as ``1e-20``.
- **Solvers apply their own tolerances.** A pass bound is met within the
  solver's feasibility tolerance and computed with its own conventions, so a
  plan can exceed the nominal bound slightly, and a near-exact bound can make
  a later pass infeasible. If a later pass is reported infeasible with a
  near-zero tolerance, use a small positive one (e.g. ``1e-6``).
- **Trade-off curves.** A tolerance turns the hierarchy into the
  epsilon-constraint method: sweeping ``.objreltol`` of the first objective
  traces the trade-off curve between it and the next one. Each point is only
  weakly efficient, since the first objective is capped, not re-optimized;
  adding a copy of it as the last (lowest-priority) objective makes the
  points efficient. See
  `AMPL Colab notebooks on Pareto fronts <https://colab.ampl.com/tags/pareto-front.html>`__.


Nonlinear multiple objectives
**************************************

In native mode (``obj:multi=1``), many solvers only accept linear objectives.
To use it for nonlinear objectives, set ``cvt:nlobj=0 cvt:quadobj=0``
which moves nonlinear terms into constraints.


.. _objective-specific-options:


Options for each objective
********************************************************

You can specify options for each objective by creating
a suffix in AMPL with the name starting with ``option_``
followed by an :ref:`option name <ampl-solver-options>`.
This should be set using the suffix notation; see the example below
where we set a time limit and a mip gap for each objective.

Objective-specific options are supported in the
:ref:`multi-objective emulator <multipleObjectives>`
(option ``obj:multi=2``), as well as natively in Gurobi
(option ``obj:multi=1``, since AMPL/Gurobi 20251211).
Set ``obj:multi:options=0`` to ignore
objective-specific option suffixes.

.. note::
   Gurobi options `obj:N:priority` do not currently work
   with objective options; use suffixes `.objpriority` as in the
   example below.

.. note::
   If the suffix value is not set, the default value
   or the initially specified value
   of the option will be used. Also note that AMPL suffixes
   have default value 0 which means 'not set'.
   Thus, to set a real-valued objective-specific
   option to 0, use a very small value, such as
   ``let _obj[1].option_mipgap := 1e-20;``. For integer-valued options
   this is currently not possible.


.. code-block:: ampl

    suffix objpriority integer; # Declare the suffixes
    suffix option_timelimit;    # Use single-word option alias
    suffix option_mipgap;


    option gurobi_options "obj:multi=2 timelimit=60 mipgap=0.00001";

    minimize total_cost {s in 1..3}:
       sum {j in FOOD} cost[s,j] * Buy[j] suffix objpriority s;
    minimize total_number:  sum {j in FOOD} Buy[j];

    # Note an alternative way to set the objective priority
    let total_number.objpriority := 3;

 
    for {s in 1..3} {
        let total_cost[s].option_timelimit := s*10;
        let total_cost[s].option_mipgap := s*0.01;
    }
    let total_number.option_timelimit := 30;

    # Note that objectives total_cost[3] and total_number will be blended
    # (since they have the same priority),  but there is no conflict
    # because the timelimit is set to 30 for both, and mipgap is set
    # only for total_cost[3], see below for explanation.



If multiple objectives have the same priority, they are blended together.
When objectives are part of the same blended group, the driver will reject configurations where 
different option values are specified for these objectives. All objectives within a blended group 
must share identical option values (or not be set).

See
`AMPL Colab notebooks with multi-objective options <https://colab.ampl.com/tags/multi-objective-options.html>`__
for examples.


Ordinary gap and bound suffixes
**************************************

For a single optimized objective, ``mip:return_gap`` requests
``.relmipgap`` / ``.absmipgap`` and ``mip:bestbound`` requests ``.bestbound``.
With multiple active objectives, use the pass statistics below instead:
ordinary suffixes do not provide independent certificates for every objective.
Some paths report a current/final-pass scalar on the first objective and leave
other entries zero. Those zeros do not establish a zero gap or a valid bound.

A pass bound applies to its weighted objective. When several objectives share
one pass, it is not generally a bound on each original objective separately.
Unavailability must not be interpreted as optimality.


.. _multiobj-pass-statistics:

Multi-objective solve statistics (experimental)
***********************************************************************************

Use the outputs according to their scope:

.. list-table::
   :header-rows: 1
   :widths: 28 42 30

   * - Output
     - Meaning
     - Availability
   * - Ordinary objective ``.relmipgap``, ``.absmipgap``, ``.bestbound``
     - Useful for a single optimized objective. In multiobjective mode,
       entries can contain a current/final-pass scalar or padding.
     - Do not infer multiobjective optimality from zero entries.
   * - Ordinary problem-level gap/bound suffixes
     - Scalar backend information; not a collection of per-pass certificates.
     - Not a substitute for multiobjective pass statistics.
   * - ``.objpass_*``
     - Statistics for the weighted optimization pass associated with an objective.
     - Emulation collects them; native availability depends on the driver.

With option `obj:multi:stats=1`, the following suffixes are available
after a multi-objective solve:

- `.objpass`: the index of the optimization pass where this objective was optimized
  (0 if not optimized due to a time limit, infeasibility, etc. in previous passes);
- `.objpass_result`: the `solve_result_num` value of the pass of this objective,
  -1 if not optimized;
- `.objpass_(runtime/mipgap/objval/objbound/work/itercount/nodecount/opennodecount)`:
  statistics of the pass of this objective, where available from the solver.

Objectives blended into the same pass receive the same values.
For example, a lower bound of 14 on ``2*f + g`` is a certificate for that
weighted sum, not separate bounds of 14 on ``f`` and ``g``. Even a zero pass
gap does not mean either component was independently minimized. Evaluate
original objectives normally after solving and use ``.objpass`` to identify
which objectives share statistics.

In Gurobi native mode (`obj:multi=1`), the statistics are as described
at `Gurobi multi-objective documentation <https://docs.gurobi.com/projects/optimizer/en/current/reference/attributes/multiobjective.html>`__.
In CPLEX native mode, they are taken from CPLEX's per-subproblem information
(`.objpass_work` is deterministic time in ticks; `.objpass_mipgap`
is computed from the objective value and the best bound,
NaN for continuous passes where CPLEX reports no bound).
Other solvers with native multi-objective support
do not provide pass statistics; use `obj:multi=2` to obtain them via emulation. Read the driver's warning and capability
information: absent output/default zero is not a reported metric. In emulation,
unavailable numeric statistics are NaN as described below.

In all modes, the pass objective is the weighted sum of the objectives
of that priority level, optimized in the sense of the first objective.
With the default `obj:multi:weight=2`, the weight of an objective
of the opposite sense is negated; with `obj:multi:weight=1`,
the weights are applied as given.

With MP's multi-objective emulation (`obj:multi=2`, or `obj:multi=1`
for solvers without native support), the statistics are collected
after each pass:

- `.objpass_objval` and `.objpass_objbound` refer to the pass objective,
  as in the native modes. For example, the value of a minimization objective
  is negated if the first objective is maximized. The original objectives
  can be evaluated as usual after the solve;
- `.objpass_runtime` is the wall-clock time of the pass, including
  the model modifications and solution retrieval between passes;
- numeric statistics not available for a pass, including passes
  not run, are reported as NaN.



Examples
**************************************

See the :ref:`multiple objectives feature <multipleObjectives>` and
`Multi-objective AMPL Colab notebooks <https://colab.ampl.com/tags/multi-objective.html>`__
for examples.
