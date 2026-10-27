
.. _numerical_accuracy:

Numerical accuracy
------------------------

Mathematical Programming solvers typically work with finite-precision numbers, which
leads to concerns on numerical stability.
While general numerical stability is a well-established topic, below
we highlight some common points.


Tight coefficients
*****************************

Try & keep the model's coefficients
(constraint matrix, objective)
close to :math:`\pm1` by rescaling your data.


.. _tight_bounds:

Importance of tight bounds: "big-M" constraints
**************************************************

Try & keep your variable and constraint bounds tight.

For example, the MP library simplifies relational operators into "indicator" constraints.
Solvers natively supporting indicators, usually handle them in a numerically stable way.
This is about numerical stability, not speed: a big-M row with a tight
constant contributes LP relaxation strength that a native indicator may not,
so an indicator formulation is not automatically faster
(see :ref:`indicator-vs-bigm`).
Otherwise, they have to be linearized by the so-called "big-M" constraints, for example,
an implication

.. code-block:: ampl

      b==0 ==> x<=0;

may be linearized as

.. code-block:: ampl

      x <= upper_bound(x) * b;

with the big-M constant taken as the upper bound on :math:`x`. Moreover,
such a linearizaton can be performed to compute a continuous relaxation
of the MIP model, facilitating faster solving.

Thus, big-M constraints require finite bounds on participating variables.
These are the **declared** bounds, propagated through the expressions;
see :ref:`declared-bounds-bigm` for how to make a bound implied by other
constraints visible to MP.
They should be as tight as possible, ideally between :math:`\pm10^4`.
In any case, for numerical stability these bounds should
not exceed the reciprocal of the integrality tolerance
(:ref:`option <solver-options>` *inttol*). A default
big-M value can be set with the option *cvt:bigM* (use with caution).

In some cases, a different formulation can help automatic detection
of bounds. In the below example, the ``==> / else`` operator
cannot easily do this for variables ``DL``:

.. code-block:: ampl

    var DL {l in Lset, s in Sset};
    subject to setDL {l in Lset, s in Sset}:
       L[l,s] > K[l,s] ==> DL[l,s] = L[l,s] else DL[l,s] = 0;

Reformulating this constraint using ``if-then`` enables automatic
bound deduction:

.. code-block:: ampl

    subject to setDL {l in Lset, s in Sset}:
       DL[l,s] = if L[l,s] > K[l,s] then L[l,s];


.. _counting-continuous-changes:

Counting changes of a continuous variable
*******************************************

A variable declared in a finite set, such as ``var A {T} in TAPS``, is still
a continuous variable for the solver, returned only up to its feasibility
tolerance. Two periods at the same value can come back as ``1.025`` and
``1.0249994``, and an exact comparison ``A[t] != A[t-1]`` then counts a change
that did not happen. Compare the move with half a step of the value ladder
instead:

.. code-block:: ampl

    # not this
    subject to Few_Moves: atmost N {t in T: t > 1} (A[t] != A[t-1]);
    # but this: real moves are at least one step, noise is ~1e-6
    subject to Few_Moves: atmost N {t in T: t > 1} (abs(A[t] - A[t-1]) >= step/2);

``!=`` is safe on variables declared ``integer``.


Piecewise-linear functions
*****************************

Piecewise-linear expressions can be modeled in AMPL directly, or arise from
approximations of other functions. Solvers which support PL expressions,
usually handle them algorithmically in a numerically stable way. Otherwise,
if PL expressions are linearized, it is recommended to have the argument
and result variables bounded in :math:`\pm10^4` (the tighter the better;
for :ref:`approximated nonlinear functions <nonlinear-functions>`,
hard bounds of up to :math:`\pm10^6` are imposed by default).
The MIP linearization of a nonconvex PL expression via SOS2 constraints
needs finite bounds on the argument: an unbounded argument is an error,
unless a limit is set with the option *cvt:pl:domain* (use with caution:
the solver's feasibility tolerance is amplified by up to this value,
which can lead to violated PL expressions). With native SOS2 support
in the solver, unbounded sides of the argument are modeled exactly
by ray variables. The stability can be improved
in some cases by decreasing integer tolerance, Gurobi's *intfocus* and
*numfocus* options, switching off presolve in the solver, and other tuning measures.

