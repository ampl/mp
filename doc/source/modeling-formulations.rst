.. _modeling-formulations:

Choosing between formulations
--------------------------------

MP lets you state logic directly (``==>``, ``or``, ``atmost``, ``count``,
set membership, ``if-then``), and every such model also has hand-written
alternatives: binaries with big-M constraints, or a reformulation from the
literature. Compare alternatives that preserve the intended model, accounting
for valid bounds, tolerances and any approximations. They can differ greatly in
solving time, and **which one is faster is often not predictable from the
model text**. This section collects rules that held up in measurement on
scheduling, energy, network design and disjunctive test problems, and cheap
checks that tell you which case you are in before a long benchmark.
For the general topic see also Vielma (2015) and Klotz & Newman (2013),
cited in :ref:`efficiency`.

Terms used below: the *root bound* is the optimum of the continuous (LP)
relaxation; a formulation is *stronger* when its root bound is closer to the
integer optimum.


.. _root-bound-check:

Compare root bounds before benchmarking
*******************************************

Solving the continuous relaxation can be a cheaper first comparison and shows
whether a reformulation improves the relaxation bound. With logical
constraints in the model, two details matter:

- Use the solver option *alg:relax* (``relax=1``). It relaxes integrality of
  **all** variables the solver receives, including the binaries MP creates.
  AMPL's ``option relax_integrality 1`` relaxes only the model's own integer
  variables, so with MP-converted logic it does not give the root bound.
- Read the objective value **the solver reports**. After a relaxed solve,
  AMPL recomputes the objective from the returned variables and evaluates
  logical expressions (``if-then``, ``==>``) exactly at the fractional point,
  which can give a value above the true optimum of a minimization. The
  :ref:`automatic solution check <solution-check>` warns about this.

.. tab-set::

   .. tab-item:: Python

      .. code-block:: python

         import re

         ampl.solve(solver="highs", highs_options="relax=1 chk:mode=0")
         message = ampl.get_value("solve_message")
         root_bound = float(re.findall(r"objective ([-+0-9.eE]+)", message)[-1])

   .. tab-item:: AMPL

      .. code-block:: ampl

         ampl: option highs_options 'relax=1 chk:mode=0';
         ampl: solve;
         HiGHS 1.15.1: optimal solution; objective 711143.2092

Use root bounds to understand a formulation, not to predict its solve time.
Equal root bounds can accompany very different search behavior
(see :ref:`tighter-not-faster`). Follow the relaxation check with fixed-budget
solves on matching data and seeds. Compare feasible objectives, remaining gaps,
time to solutions where available, and limits/failures as well as solved-case
times. Keep threads, limits, tolerances and memory controls fixed, and check
returned points against the original model.


.. _indicator-vs-bigm:

Indicator constraints or big-M
*********************************

Start by expressing the implication directly in AMPL, for example
``b = 1 ==> constraint``. MP passes it to a solver's native indicator handler
where supported and selected, or reformulates it otherwise. This keeps one
semantic model while allowing a controlled comparison of treatments.

Solve first with the driver's defaults. Then solve the same model and data
with the following MP options to request linearization of linear indicators:

.. code-block:: text

    acc:indeq=0 acc:indle=0 acc:indge=0

Place these controls in the MP driver's option string; with an embedded MP2NL
bridge use ``mp_options`` (see :ref:`mp2nl-option-routing`). Use the same time,
thread and numerical settings in both runs. Check the generated constraint
counts or exported model to confirm which treatment was used. See
:ref:`acceptance-options`; ``acc:_all=0`` is a broader experiment, changing other
constructs too. These options do not linearize every kind of logical expression
or request quadratic-indicator linearization.

Provide tight, justified variable bounds: MP sizes big-M constants from bounds
(see :ref:`declared-bounds-bigm`). Native indicators avoid requiring a manually
chosen M, but they are not universally faster. A valid tight big-M formulation
can strengthen the relaxation (Bonami, Lodi, Tramontani & Wiese 2015); native
handling and solver presolve can lead to different search behavior. Bounds can
matter to both treatments, so measure rather than concluding that tightening
bounds cannot help indicators.

Write a separate formulation when it adds structural knowledge that the
option-based comparison does not: a tighter data-derived M, convex-hull
inequalities, perspective constraints or a scheduling formulation from the
literature. A hand-written big-M need not equal MP's bound-derived one. Check
mathematical equivalence, variable counts and original-model residuals before
attributing a timing difference to indicators alone.


Logic in implications, and facets
***********************************

AMPL transmits ``b = 1 ==> C`` as a disjunction, and MP recovers indicator
constraints from it, including for a conjunction of comparisons
(``b = 1 ==> forall {s in S} x[s] = 1``), see option *cvt:pre:disj2ind* in
:ref:`conditional-operators`. On
a unit commitment test, a minimum up-time written with a ``forall``
consequent gave exactly the same root bound as one implication per period.

Where the literature has a stronger inequality for the same structure, it is
worth writing it even in an otherwise logical model. For minimum up/down
times the window inequalities of Rajan & Takriti (2005) are facet-defining:

.. code-block:: ampl

    # implication form
    subject to Min_Up {u in U, t in T}:
       Start[u,t] = 1 ==> forall {s in t..t+minup[u]-1: s in T} On[u,s] = 1;

    # facet-defining form (Rajan & Takriti 2005)
    subject to Min_Up_RT {u in U, t in T}:
       sum {s in t-minup[u]+1..t: s in T} Start[u,s] <= On[u,t];

On one hydro scheduling model, the best formulation was a *hybrid*: MP's
set-membership declaration for the forbidden operating zones and the facet
inequalities for the dwell times. It beat both the purely logical and the
purely hand-written formulation at every size. Restate what MP states well,
and keep what the literature has proven strongest.


Known formulations for common structures
*******************************************

Before writing a mechanical big-M translation, check whether a better
formulation exists for the structure:

.. list-table::
   :header-rows: 1
   :widths: 25 30 25 20

   * - Structure
     - Formulation
     - Reference
     - When it does not pay
   * - scheduling with start times and resources
     - time-indexed: one binary per (job, start period)
     - Pritsker, Watters & Wolfe 1969; Sousa & Wolsey 1992
     - long horizons with fine time steps
   * - ``alldiff`` over a small domain
     - assignment: one binary per (item, value)
     - Vielma 2015
     - domains much larger than the number of items
   * - semicontinuous variable with a convex cost
     - perspective reformulation (:ref:`perspective`)
     - Frangioni & Gentile 2006; Günlük & Linderoth 2010
     - linear or concave costs
   * - product of a discrete choice and a continuous variable
     - disaggregation per choice (convex hull of binary times bounded continuous)
     - Glover 1975; Raghunathan 2013
     - many choices times many products
   * - minimum up/down times
     - facet-defining window inequalities
     - Rajan & Takriti 2005
     - the constraint does not bind
   * - disjunctions in general
     - convex hull or big-M
     - Balas 1985; Grossmann & Trespalacios 2013
     - see :ref:`tighter-not-faster`

Count the variables and constraints a formulation needs at your instance
sizes before building it; one candidate in our tests needed 13 to 22 times
more binaries than the formulation it was meant to replace.


.. _perspective:

Perspective reformulation
***************************

For a convex cost :math:`f(x)` with :math:`f(0)=0`, paid only while a unit is
on (:math:`z=1`), the textbook rows :math:`\ell z \le x \le u z`,
:math:`\mathit{cost} \ge f(x)` are exact but have a weak relaxation. The
perspective :math:`\mathit{cost} \ge z\,f(x/z)` together with the on/off
bounds describes the convex hull of "off, or on and paying :math:`f`"
(Frangioni & Gentile 2006). It is rarely written literally:

- **Piecewise-linear convex** :math:`f = \max_k (s_k x + c_k)`, often written
  with AMPL's ``<<...>>`` notation (:ref:`piecewise_linear_modeling`): multiply
  every **intercept** by :math:`z`. The rows stay linear:

  .. code-block:: ampl

      param c {k in 1..nseg} := sum {m in 1..k-1} (s[m] - s[m+1]) * brk[m];  # c[1] = 0
      subject to Cost_Persp {k in 1..nseg}: Cost >= s[k] * X + c[k] * Z;

- **Quadratic** :math:`f = q x^2`: a rotated second-order cone,
  ``q * X^2 <= Cost * Z`` (see :ref:`socp`).
- **With MP's semicontinuous declaration** ``var X in {0} union interval[lo, hi]``
  there is no :math:`z` in the model. Gate the intercepts with a threshold
  inside the gap :math:`(0, lo)`; MP decides the comparison with the binary it
  creates for the declaration (option *cvt:pre:vbnd*, see
  :ref:`conditional-operators`), so the result is the same perspective:

  .. code-block:: ampl

      subject to Cost_Seg1: Cost >= s[1] * X;
      subject to Cost_Persp {k in 2..nseg}: Cost >= s[k] * X + (if X >= lo/2 then c[k]);

On a unit commitment test the perspective closed 99% of the root gap of the
textbook formulation. Whether that makes the solve faster depends on the
solver; on one solver the plain piecewise-linear form was still the fastest.
Do not convert perspective rows into indicators (``Z = 1 ==> Cost >= ...``):
that discards the strengthening at fractional :math:`z`. A fixed cost
:math:`\mathit{fix} \cdot z` already is a perspective. See the notebooks
tagged `perspective <https://colab.ampl.com/tags/perspective.html>`__.


Set membership or disjunction
*******************************

Writing ``var X in SOMESET`` or a disjunction ``or {k in SOMESET} X = k``
is a formulation choice with no known sign. On two models of the same
shape, a product of a discrete-valued variable and a continuous one, the
two forms won in opposite directions, by a factor of 4 to 600 in the
remaining gap. Until the set membership rework is complete, see the
recommendation in :ref:`set-membership`, and measure both forms on your
model.


.. _tighter-not-faster:

A tighter formulation is not always faster
*********************************************

Strength trades against size, and the balance shifts with scale and with the
solver:

- Convex hull against big-M crossed over with instance size on the same
  disjunctive test problem.
- On two catalogue-selection models, big-M and the convex hull had
  identical root bounds; the hull won on one and lost all 12 comparisons
  (up to 15 times slower) on the other. The root bound does not predict the
  winner.
- An exact tightening of a variable window helped two solvers and made a
  third one stall on a worse incumbent, because it removed the slack the
  solver's primal heuristics were using.

So write both formulations and benchmark them, after making the comparison
fair: size the big-M constant honestly (a hull that beats a deliberately
loose big-M proves nothing), and check that the formulation really is what
you intended (count variables and rows, compare root bounds); the checklist
is in :ref:`before-benchmarking`.


.. _substitute-or-name:

Nonlinear models: substitute or name
***************************************

The same choice exists for nonlinear models. A quantity defined by an
equality can be **substituted** (AMPL defined variables ``var x = expr``, or
``option substout 1``) or kept as a **named variable** with its defining
constraint.

- Substitution removes the variable and the row. On a large power network
  NLP it gave 25% fewer variables and rows and a 1.7 times faster solve,
  with identical optima.
- Keep the named variable when the name lets MP recognize structure (a cone
  ``Wc^2 + Ws^2 <= Wi * Wj`` over lifted variables is recognized; with the
  lifts inlined it is a degree-4 expression), when the quantity needs a bound
  that its expression does not imply (defined variables carry no bounds; see
  the example in :ref:`subexpression-bounds`), or when the defining
  constraint's dual is needed.
- Coordinates are a choice too: a polynomial form (for example rectangular
  voltages) suits solvers that only accept quadratic nonconvexity and conic
  relaxations; the original transcendental form has fewer variables and can
  be better for solvers that handle general nonlinear functions natively.
  Published comparisons disagree, so benchmark the alternatives on the
  solvers you use.
- A constraint that is an identity of the model's own definitions can still
  give a global solver the convex bound it lacks, see
  :ref:`redundant-convex-constraints`.

Keep both forms behind a parameter so the choice stays testable:

.. code-block:: ampl

    param defvars binary default 1;
    var FlowF {K} >= -cap[k], <= cap[k];
    var Flow {k in K} = if defvars = 0 then FlowF[k] else <expression>;
    subject to Flow_Def {k in K: defvars = 0}: FlowF[k] = <expression>;


Comparing formulations fairly
*******************************

Validate a reformulation against the original on identical data and more
than one solver before trusting its timings, and compare times only over
instances both formulations solved. See :ref:`reading-solver-results` and
:ref:`comparing-formulations`.


References
************

- Balas, E. (1985). Disjunctive programming and a hierarchy of relaxations for
  discrete optimization problems. *SIAM J. Algebraic Discrete Methods* 6(3), 466-486.
- Bonami, P., Lodi, A., Tramontani, A., Wiese, S. (2015). On mathematical
  programming with indicator constraints. *Math. Programming* 151(1), 191-223.
- Frangioni, A., Gentile, C. (2006). Perspective cuts for a class of convex
  0-1 mixed integer programs. *Math. Programming* 106(2), 225-236.
- Glover, F. (1975). Improved linear integer programming formulations of
  nonlinear integer problems. *Management Science* 22(4), 455-460.
- Grossmann, I.E., Trespalacios, F. (2013). Systematic modeling of
  discrete-continuous optimization models through generalized disjunctive
  programming. *AIChE J.* 59(9), 3276-3295.
- Günlük, O., Linderoth, J. (2010). Perspective reformulations of mixed
  integer nonlinear programs with indicator variables. *Math. Programming*
  124, 183-205.
- Pritsker, A.A.B., Watters, L.J., Wolfe, P.M. (1969). Multiproject scheduling
  with limited resources: a zero-one programming approach. *Management
  Science* 16(1), 93-108.
- Raghunathan, A.U. (2013). Global optimization of nonlinear network design.
  *SIAM J. Optim.* 23(1), 268-295.
- Rajan, D., Takriti, S. (2005). Minimum up/down polytopes of the unit
  commitment problem with start-up costs. IBM Research Report RC23628.
- Sousa, J.P., Wolsey, L.A. (1992). A time indexed formulation of
  non-preemptive single machine scheduling problems. *Math. Programming* 54,
  353-367.
- Vielma, J.P. (2015). Mixed integer linear programming formulation
  techniques. *SIAM Review* 57(1), 3-57.

Worked examples: the notebooks tagged
`formulation-comparison <https://colab.ampl.com/tags/formulation-comparison.html>`__.
