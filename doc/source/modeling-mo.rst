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

Objectives can be blended, or hierarchical (Lexicographical).


Blended objectives
************************************************

By default, all objectives are blended together
(summed up with the corresponding signs.)
Suffixes ``.objweight`` can be used to change the individual weights
and objective senses, according to the option ``obj:multi:weight``.


Lexicographical objectives
********************************************************

To apply hierarchical optimization, use suffix ``.objpriority``,
as described in the ``obj:multi`` option description.

.. code-block:: ampl

    suffix objpriority;

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
(option ``obj:multi=1``). Set ``obj:multi:options=0`` to ignore
objective-specific option suffixes.

Note that if the suffix value is not set, the default value
or the initially specified value
of the option will be used. Also note that AMPL suffixes
have default value 0 which means 'not set'.
Thus, to set a real-valued objective-specific
option to 0, use a very small value, such
``let _obj[1].option_mipgap := 1e-20;``. For integer-valued options
this is currently not possible.


.. code-block:: ampl

    suffix objpriority;
    suffix option_timelimit;    # Use single-word option alias
    suffix option_mipgap;


    option gurobi_options "timelimit=60 mipgap=0.00001";

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

    # note that objectives total_cost[3] and total_number will be blended,  but there
    # is no conflict because the timelimit is set to 30 for both, and mipgap is set
    # only for total_cost[3].



If multiple objectives have the same priority, they are blended together.
When objectives are part of the same blended group, the driver will reject configurations where 
different option values are specified for these objectives. All objectives within a blended group 
must share identical option values or at most one should have the value specified.

See
`AMPL Colab notebooks with multi-objective options <https://colab.ampl.com/tags/multi-objective-options.html>`__
for examples.


Examples
**************************************

See the :ref:`important-features` and
`Multi-objective AMPL Colab notebooks <https://colab.ampl.com/tags/multi-objective.html>`__
for examples.
