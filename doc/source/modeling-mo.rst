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
mode (``obj:multi=1``), in particular for LPs. Consult the solver documentation.

Examples
**************************************

See the :ref:`important-features` and
`Multi-objective AMPL Colab notebooks <https://colab.ampl.com/tags/multiple-objectives.html>`_
for examples.
