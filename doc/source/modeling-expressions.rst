
.. _expressions_supported:

Expressions supported
---------------------

The MP solver interface library works with existing AMPL syntaxes, but allows them to be used in more general ways, or with a greater variety of solvers.

In many cases, an extension results from allowing variables to appear in more general contexts, such as with conditional, logical, or counting operators. Other extensions are enabled by providing more powerful transformations, particularly to linear or quadratic equivalents, and by providing support for extensions that are native to some solvers. A few extensions are already handled in the AMPL language translator, and are included here for completeness.

In the syntax summaries below, there are two main kinds of entities, representing *numerical expressions* and *constraints:*

- **expr**
     represents any expression that evaluates to a number. Unless otherwise indicated, it may contain variables. It may be built from familiar arithmetic operators, but also from other operators or functions that return numerical values.

- **constr**
     represents a constraint of the model, which may evaluate to true or false
     depending on the values of variables that it contains. It may be built from the
     familiar relational operators ``>(=)``, ``<(=)``, and ``=``, but also from other
     operators such as ``or`` and ``alldiff`` that create constraints.

The return value of an operator or function is also one of the above, as indicated by *expr-valued* or *constr-valued* at the beginning of each syntax summary. Thus it is possible to build up complex combinations of operators and functions of various kinds; for example,

.. code-block:: ampl

        x<=0 or (y!=2  ==>
                 x<=-5 else
                       max((x+1)*x*z, y, y-z)<=3 and exp(y)<=12);

AMPL represents these combinations as expression trees, which are sent to MP-based solver interfaces to be processed as solvers require.

Indexing over sets is a common feature of AMPL expressions. The examples below use two kinds of indexing expressions, which are represented in the syntax summaries as follows:

- { indexing }
    This is the regular sort of AMPL indexing expression, as used in defining numerous AMPL entities such as parameters, variables, constraints, and summations. It is described in the `AMPL book <https://ampl.com/resources/the-ampl-book/>`_ beginning with `Section 5.5 Indexing expressions <https://ampl.com/BOOK/CHAPTERS/08-sets1.pdf#page=7>`_ and continuing with `Chapter 6. Compound Sets and Indexing <https://ampl.com/BOOK/CHAPTERS/09-sets2.pdf>`_. Followed by an *expr* or *constr*, an indexing expression specifies a list of expressions or constraints to which an operator applies; for example,
    ::

        max {n in NODE} weight[t,n] * Use[n]
        forall {p in PROD} Trans[i,j,p] = 0

- ( expr-list )
    This is a parenthesized, comma-separated list of entries that represent numerical values. Each entry may have the form *expr* or *{indexing} expr*, or recursively *{indexing} ( expr-list )*. For example,
    ::

        max (cost["BRO"],cost["CAU"],cost["BRU"])
        max ({f in FOOD} cost[f], 10.0)
        max ({n in NUTR} (lim_nutr[n], {f in FOOD} amt[n,f]))

As seen in the case of ``max`` above, certain operators can be used with either the ``{indexing} expr`` or the ``(expr-list)`` form.

Due to the generality of the operators recognized by the MP interface, it is possible to express constraints that do not define a closed feasible region. For example,

.. code-block:: ampl

        x > 5
        not (x >= 5 and y >= 15)
        x = 0 ==> z = 0 else z = 1

An optimum is not guaranteed to exist over a non-closed region.
Thus where necessary, the MP interface constructs an approximate closed region by
use of a small tolerance. For example, if x is minimized subject to x > 5, then any
x value greater than 5 is not minimal, and any x value less than or equal to 5 is
not feasible. Thus, to insure that a minimum is well defined, the constraint must
be changed to x >= 5 + eps for some small constant eps. Each solver has its own
default value of the eps constant, which can be adjusted through
an :ref:`option setting <solver-options>`.


Conditional operators
***********************************

- if *constr* then *expr1* [else *expr2*]
    *expr-valued:* When *constr* is true, takes the value of *expr1*.
    When *constr* is false, takes the value of *expr2*, or 0 if the ``else`` phrase is omitted.

In the special case where there are no variables in the *constr*, the value of this expression
can be determined as either *expr1* or *expr2* (or 0) before the problem is sent to the solver.
But in general, the value of expression depends upon how the solver assigns values to the
variables in the *constr*, and so AMPL must send the entire expression to the solver interface for processing.

.. code-block:: ampl

       minimize TotalCost:
          sum {j in JOBS, k in MACHINES}
             if MachineForJob[j] = k then cost[j,k];

.. code-block:: ampl

       subject to Balance {p in PROD, t in TIME}:
          Make[p,t] + (if t = 0 then inv0[p] else Inv[p,t-1])
             = Sell[p,t] + Inv[p,t];

- *constr1* ==> *constr2* [else *constr3*]
    *constr-valued:* Satistifed when *constr1* is true and *constr2* is true,
    or when *constr1* is false [and also *constr3* is true, if present].
- *constr2* <== *constr1*
    *constr-valued:* Satistifed when *constr1* is true and *constr2* is true,
    or when *constr1* is false.
- *constr1* <==> *constr2*
    *constr-valued:* Satisfied if *constr1* and *constr2* are both true or both false.

The conditional expression *constr1* ==> *constr2* can be thought of as saying that
*constr1* implies *constr2*, or equivalently that if *constr1* then *constr2*. In the
special case where *constr1* is of the form *binary-var* = 0 or *binary-var* = 1, these
are "indicator" constraints that can be handled natively by some solvers. Otherwise,
they are transformed to simpler constraints that use relational operators. The other
cases are treated similarly.

.. code-block:: ampl

    subject to Multi_Min_Ship {i in ORIG, j in DEST}:
       sum {p in PROD} Trans[i,j,p] > 0 ==>
          minload <= sum {p in PROD} Trans[i,j,p] <= limit[i,j];

.. code-block:: ampl

    subject to Least_Use {j in SCHEDS}:
       Use[j] = 1 ==> Work[j] >= least_assign else Work[j] = 0;


Logical operators
***********************************

- *constr1* or *constr2*
    *constr-valued:* Satisfied when *constr1* is true or *constr2* is true.
- *constr1* and *constr2*
    *constr-valued:* Satisfied when *constr1* is true and *constr2* is true.
- not *constr*
    *constr-valued:* Satisfied when *constr* is false.

Expressions using these operators are transformed to use Gurobi's native AND
and OR "general constraints" when possible. In other cases, they are
transformed to simpler constraints that use relational operators.

.. code-block:: ampl

    subj to NoPersonIsolated
             {l in TYPES['loc'], r in TYPES['rank'], j in 1..numberGrps}:
       sum {i in LOCRANK[l,r]} Assign[i,j] = 0 or
       sum {i in LOCRANK[l,r]} Assign[i,j] +
          sum {a in ADJACENT[r]} sum {i in LOCRANK[l,a]} Assign[i,j] >= 2;

.. code-block:: ampl

    subj to No_Conflict {i1 in JOBS, i2 in JOBS: ord(i1) < ord(i2)}:
       Start[i2] >= Start[i1] + t_offset[i1,i2]  or
       Start[i1] >= Start[i2] + t_offset[i2,i1];

.. code-block:: ampl

    subject to Least_Use {j in SCHEDS}:
       Work[j] = 0 or least_assign <= Work[j] <= max {i in SHIFT_LIST[j]} required[i];

.. code-block:: ampl

    subj to EntRem {t in 1..numTanks}:
       Entry[t] + minTime[t] <= Removal[t] and
       Entry[t] + maxTime[t] >= Removal[t];

- exists {indexing} *constr*
    *constr-valued:* Satisfied when at least one of the *constr* operands is true.
- forall {indexing} *constr*
    *constr-valued:* Satisfied when all of the *constr* operands are true.

The ``exists`` and ``forall`` operators are the iterated forms of ``or`` and ``and``, respectively.

.. code-block:: ampl

    minimize Total_Cost:
       sum {p in PRODUCTS, (i,j) in ARCS} var_cost[p,i,j] * Flow[p,i,j] +
       sum {(i,j) in ARCS} if exists {p in PRODUCTS} Flow[p,i,j] > 0 then fix_cost[i,j];

.. code-block:: ampl

    subject to Multi {i in ORIG, j in DEST}:
       forall {p in PROD} Trans[i,j,p] = 0  or
       minload <= sum {p in PROD} Trans[i,j,p] <= limit[i,j];

.. code-block:: ampl

    subj to HostNever {j in BOATS}:
       isH[j] = 1 ==> forall {t in TIMES} H[j,t] = j;


.. _piecewise_linear_modeling:

Piecewise-linear expressions
***********************************

- abs (*expr*)
    *expr-valued:* Equals *expr* when ≥ 0, or *-expr* when < 0.
- min {indexing} *expr*
    *expr-valued:* Equals the smallest value among the *expr* operands.
- min ( expr-list )
    *expr-valued:* Equals the smallest value among all of the operands in the *expr-list*.
- max {indexing} *expr*
    *expr-valued:* Equals the largest value among the *expr* operands.
- max ( expr-list )
    *expr-valued:* Equals the largest value among all of the operands in the *expr-list*.

Expressions using these operators are transformed to use Gurobi's native ABS, MIN, and MAX "general constraints" when possible. In other cases, they are transformed to simpler constraints that use relational operators, and in particular are linearized where all of the operands are linear.

.. code-block:: ampl

    maximize Total_Profit:
       sum {p in PROD, t in 1..T} revenue[p,t]*Sell[p,t] -
       sum {t in 1..T} time_penalty[t] * abs(Use[t] - avail_min[t]);

.. code-block:: ampl

    minimize Max_Cost:
       max {i in PEOPLE} sum {j in PROJECTS} cost[i,j] * Assign[i,j];

.. code-block:: ampl

    maximize WeightSum:
       sum {t in TRAJ} max {n in NODE} weight[t,n] * Use[n];

- << *slope-list*; *breakpoint-list* >> var
    *expr-valued:* Computes a piecewise-linear function of a single variable; see
    `Chapter 17. Piecewise-Linear Programs <https://ampl.com/BOOK/CHAPTERS/20-piecewise.pdf>`_ in
    the `AMPL book <https://ampl.com/resources/the-ampl-book/>`_ for a complete description of the
    forms that AMPL recognizes.

This piecewise-linear expression is defined by lists of ``n`` *breakpoints* and ``n+1``
*slopes*. The *var* must be a reference to a single variable.

When AMPL's option ``pl_linearize`` is at its default value of 1, AMPL linearizes these
piecewise-linear expressions, and sends the linearized versions to the solver. The linearization
is continuous where possible, in certain convex and concave cases (where the slopes are
increasing and decreasing, respectively); but in general, the linearization includes both
continuous and binary variables.

When ``pl_linearize`` is set to 0, piecewise-linear expressions are represented to the solver
in the form of expression trees. The MP-based interface transforms them to use a solver's native
methods for piecewise-linear functions (Gurobi, COPT), and linearizes them for other solvers (HiGHS).

When a piecewise-linear function is linearized (rather than being handled natively by the solver),
numerical accuracy becomes a concern. To promote numerical stability, it is recommended that
the argument and result variables be explicitly bounded within [-1e+4,+1e-4]. See more in the section
on :ref:`numerical_accuracy`.


.. code-block:: ampl

    maximize Total_Profit:
       sum {p in PROD, t in 1..T} (revenue[p,t]*Sell[p,t] -
          prodcost[p]*Make[p,t] - <<0; -backcost[p],invcost[p]>> Inv[p,t]) -
       sum {t in 1..T} <<avail_min[t]; 0,time_penalty[t]>> Use[t]
       sum {p in PROD, t in 1..T}
          <<commit[p,t]; -100000,0>> (Sell[p,t],commit[p,t]);

.. code-block:: ampl

    minimize Total_Cost:
       sum {i in ORIG, j in DEST}
          <<{p in 1..npiece[i,j]-1} limit[i,j,p];
            {p in 1..npiece[i,j]} rate[i,j,p]>> Trans[i,j];


Counting operators
***********************************

- count {indexing} *constr*
    *expr-valued:* The number of members of the indexing set such that the *constr* is satisfied.

AMPL’s ``count`` operator examines an indexed collection of constraints, and returns the number of those constraints that are satisfied. The AMPL translator instantiates the specified constraint for each member of the indexing set, and communicates all of the instantiated constraints to the solver interface; then the solver interface transforms the counting operation to a form that the solver can accept.

.. code-block:: ampl

    subject to Min_Serve {i in ORIG}:
        count {j in DEST} (Ship[i,j] >= minload) >= minserve;

- atleast k {indexing} *constr*
    *constr-valued:* Satisfied when the *constr* is satisfied for at least ``k`` members of the indexing set.
- atmost k {indexing} *constr*
    *constr-valued:* Satisfied when the *constr* is satisfied for at most ``k`` members of the indexing set.
- exactly k {indexing} *constr*
    *constr-valued:* Satisfied when the *constr* is satisfied for exactly ``k`` members of the indexing set.

``k`` must be a constant arithmetic expression that evaluates to a nonnegative integer. These operators provide easier-to-read alternatives for special cases of constraints that rely on ``count``. Compare for example the ``Min_Serve`` constraint below to the one given previously using ``count``.

.. code-block:: ampl

    subject to Min_Serve {i in ORIG}:
        atleast minserve {j in DEST} (Ship[i,j] >= minload);

.. code-block:: ampl

    subj to CapacityOfMachine {k in MACHINES}:
        atmost cap[k] {j in JOBS} (MachineForJob[j] = k);

- numberof *expr* in ( *expr-list* )
    *expr-valued:* The number of items in the *expr-list* having the same value as *expr*.

This operator can provide an easier-to-read alternative for a special case of count.
Compare for example the ``CapacityOfMachine`` constraint below to the one given previously
using ``atmost``.

.. code-block:: ampl

    subj to CapacityOfMachine {k in MACHINES}:
        numberof k in ({j in JOBS} MachineForJob[j]) <= cap[k];

.. code-block:: ampl

    subj to MinInGrpDefn {j in 1..numberGrps}:
       MinInGrp <= numberof j in ({i in PEOPLE} Assign[i]);


Relational and comparison operators
***********************************

- expr1 > expr2, expr1 >= expr2
    *constr-valued:* Satisfied when *expr1* is strictly greater (or equal) than *expr2*.
- expr1 < expr2, expr1 <= expr2
    *constr-valued:* Satisfied when *expr1* is strictly less (or equal) than *expr2*.
- expr1 == expr2, expr1 != expr2
    *constr-valued:* Satisfied when *expr1* does (not) equal *expr2*.

Where possible, the MP interface transforms strict operations to ones involving ``>=`` and ``<=``,
so that optimization solvers can handle them. For example, this can be done when *expr1* and
*expr2* are integer-valued, or when an expression like ``if Flow[i,j] > 0 then fixed[i,j]``
expresses a fixed cost in an objective to be minimized. Where this is not possible, a small
tolerance is introduced, as discussed in :ref:`expressions_supported`. Relational operators
require careful modeling in regard to :ref:`numerical_accuracy`.


.. code-block:: ampl

    minimize TotalCost:
       sum {p in PRODUCTS, (i,j) in ARCS} var_cost[p,i,j] * Flow[p,i,j] +
       sum {(i,j) in ARCS}
          if sum {p in PRODUCTS} Flow[p,i,j] > 0 then fix_cost[i,j];

.. code-block:: ampl

    subject to Different_Colors {(c1,c2) in Neighbors}:
       Color[c1] != Color[c2];

- alldiff {indexing} *expr*
    *constr-valued:* Satisfied when *expr* takes a different value for every member of the indexing set.

- alldiff ( expr-list )
    *constr-valued:* Satisfied when all of the items in the *expr-list* take different values.

This operator provides a much more concise alternative to specifying ``!=`` between all pairs
in a specified collection of expressions. Currently none of the MP-based solvers support this
operator natively, so the interface transforms it to a representation in terms of simpler
constraints.

.. code-block:: ampl

    subject to OnePersonPerPosition:
       alldiff {i in 1..nPeople} Pos[i];

.. code-block:: ampl

    subject to Regions {I in 1..9 by 3, J in 1..9 by 3}:
       alldiff {i in I..I+2, j in J..J+2} X[i,j];


Complementarity operator
***********************************

- *constr1* complements *constr2*
    *constr-valued:* Satisfied when both *const1* and *constr2* are satisfied, and at least one of them holds with equality. Each of *constr1* and *constr2* must have the form *expr1 <= expr2* or *expr1 >= expr2* (and the trivial special case *expr1 = expr2* is also recognized).
- *expr* complements *constr*,  *constr* complements *expr*
     *constr-valued:* Satisfied when *constr* is satisfied, and when also if *expr* is positive then *constr* holds with equality at its lower bound, or if *expr* is negative then *constr* holds with equality at its upper bound. The *constr*  must have the form *lb <= expr <= ub* or *ub >= expr >= lb* where *lb* and *ub* are lower and upper bound expressions not involving variables.

The ``complements`` operator provides a convenient, streamlined way of expressing a common kind of relationship between two single-inequality constraints, or between an expression and a double-inequality constraint. This relationship appears in the complementary slackness conditions necessary for optimality of certain optimization problems, and in equilibrium conditions for games and for various physical systems. See `Chapter 19. Complementarity Problems <https://ampl.com/BOOK/CHAPTERS/22-complement.pdf>`_ in the `AMPL book <https://ampl.com/resources/the-ampl-book/>`_ for a detailed presentation.

Certain nonlinear solvers, notably Knitro, handle complementarity constraints natively. For MP-based solvers, the interface converts uses of ``complements`` to equivalent constraints using logical operators.

.. code-block:: ampl

    subject to Pri_Compl {i in PROD}:
       Price[i] >= 0 complements
          sum {j in ACT} io[i,j] * Level[j] >= demzero[i] - demrate[i] * Price[i];

.. code-block:: ampl

    subject to Lev_Compl {j in ACT}:
       level_min[j] <= Level[j] <= level_max[j] complements
          cost[j] - sum {i in PROD} Price[i] * io[i,j];



Set membership operator
**********************************

- var *var-name* in *set-expr* ;
    Defines a variable that must be a member of a specified AMPL set,
    as given by the expression *set-expr*. All members of the set must be numbers.

This is the simplest use of ``in`` to restrict the domain of a set; more
generally, the *in set-expr* phrase may appear in any ``var`` definition
that does not contain an *=* phrase.

Before sending a problem to the solver interface, AMPL converts variable
definitions of this kind to alternative definitions that do not use the
``in`` operator. This may involve the definition of auxiliary binary
variables and additional constraints. In the usual case where *set-expr*
is a finite set, AMPL also defines suffixes ``.sos`` and ``.sosref`` which
can be used by the solver interface to recognize variables and constraints
that have been created to implement an ``in`` operator, and to support
solvers that handle arbitrary variable domains by means of
"special ordered sets of type 1". It is also possible to specify sets
that contain continuous intervals -- and hence are infinite -- by using
the AMPL expression *interval[expr1,expr2]*.

.. code-block:: ampl

    var Buy {f in FOODS} in {0,10,30,45,55};

.. code-block:: ampl

    var Ship {(i,j) in ARCS}
       in {0} union interval[min_ship,capacity[i,j]];



Quadratic and power operators
**********************************

- *expr1* * *expr2*
    *expr-valued:* Multiplication of *expr1* and *expr2*.
- *expr1* / *expr2*
    *expr-valued:* Division of *expr1* by *expr2*.
- *expr1* ^ *expr2*
    *expr-valued:* *expr1* raised to the *expr2* power, for the special cases where
    either *expr1* or *expr2* is a constant. For *expr2* positive integer, the operator
    is decomposed into quadratic constraints if the solver supports them,
    otherwise passed to the solver natively or approximated by a piecewise-linear function.

For quadratic expressions of the form *linear \* linear* and *linear^2*, the operands
are multiplied out so that coefficients of individual quadratic terms can be extracted.
If the solver natively handles quadratic terms, then the quadratic coefficients are
passed to the solver, which decides whether and how to handle them. Otherwise, quadratic
terms are linearized where possible, such as where one of the operands is a binary variable,
or approximated.

Piecewise linearization allows handling of nonconvex QP and nonlinear models
by convex MIP solvers.
For convex MIQP solvers,
to apply linearization of quadratic expressions (it is the default for linear solvers only),
use options ``cvt:quadobj=0``, ``cvt:quadcon=0``.

Other expressions involving these operators are converted, where possible, to simpler
quadratic expressions and equality constraints through the use of auxiliary variables;
then the resulting quadratic expressions and equality constraints are handled in ways
previously described. For example:

- ``(x-1)^3`` is converted to ``(x-1) * y`` with the added constraint ``y = (x-1)^2``.
- ``x * max {j in 1..n} y[j]`` is converted to ``x * z`` with the added constraint
  ``z = max {j in 1..n} y[j]``.
- ``x / sum {j in 1..n} y[j]`` is converted to ``z`` with the added constraints
  ``z * t = x``, ``t = sum {j in 1..n} y[j]``, and ``t != 0``.

.. code-block:: ampl

    subj to Eq {i in J} :
       x[i+neq] / (b[i+neq] * sum {j in J} x[j+neq] / b[j+neq]) =
          c[i] * x[i] / (40 * b[i] * sum {j in J} x[j] / b[j]);


Conic optimization
**************************************

Some solvers can handle conic constraints with tailored algorithms:
Mosek, Gurobi, COPT. Note that general non-linear solvers accept them too,
but might not provide any specialized methods.
See `conic examples <https://colab.ampl.com/tags/conic.html>`_
at Google Colab.

Second-order cone programming (SOCP)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SOCP constraints (quadratic cones) are recognized by major commercial MIP solvers
from their algebraic representations. For some representations and solvers,
MP library provides additional conversion into solver-specific conic forms. Examples:

- Standard SOC:

  .. code-block:: ampl

     x[0] >= sqrt(x[1]^2 + ... + x[n]^2);
     -5*x[0]^2 <= -x[1]^2 - ... - x[n]^2,     where x[0] >= 0;
     0.04 >= abs(x[1]);


- Rotated SOC  (where x[0], x[1] >= 0):

  .. code-block:: ampl

     3*sqrt(5*x[0]*x[1]) >= 15*sqrt(10*x[2]^2 + ... + 80*x[n]^2);
     2*x[0]*x[1] >= x[2]^2 + ... + x[n]^2;


*Note:* Mosek cannot mix SOCP and general quadratic constraints,
complaining::

  MOSEK 10.0.43:  Error type 3, MSK_RES_ERR_MIXED_CONIC_AND_NL(1501):
    The problem contains both conic and nonlinear constraints.


In this case, setting :ref:`option <solver-options>` ``cvt:socp=0`` results in second-order conic
constraints being passed to the solver as quadratics, even if
the solver has native SOCP API. This gives the solver a chance to recognize
SOCP forms in preprocessing.

Exponential cones
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Mosek natively handles exponential conic constraints.
They are recognized by MP and passed to Mosek via its API.
Example:

.. code-block:: ampl

   var q1;
   var q2 >= 0;
   var w >= 0;

   # conic constraints
   s.t. T1: 1 + b*w >= exp( q1 );
   s.t. T2: -1 + w +10*q2   <= -5 * q2 * exp( q1 / (q2*5) );

*Note:* Mosek cannot mix conic and general quadratic constraints.
Use :ref:`option <solver-options>` ``cvt:expcones=0`` to handle exponential conic
constraints in other ways, even if
the solver has native exponential conic API. For Mosek, they would be approximated
by :ref:`piecewise_linear_modeling`.

Exploring the final model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To explore the model received by the solver,
e.g., to see which constraints are received as conic vs quadratic, export the model
in one of the solver's general formats::

  ampl: option mosek_auxfiles rc;   ## To use var/con names
  ampl: option mosek_options 'writeprob=/tmp/ell.jtask'; solve;




General nonlinear functions
**********************************

Supported functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- log (*expr*), log10 (*expr*)
    *expr-valued:* The natural and base-10 logarithms of *expr*.
- exp (*expr*)
    *expr-valued:* The base of the natural logarithm (e) raised to the power *expr*.
- sin (*expr*), cos (*expr*), tan (*expr*), asin (*expr*), acos (*expr*), atan (*expr*)
    *expr-valued:* The sine, cosine, tangent of *expr* and the corresponding inverse functions.
- sinh (*expr*), cosh (*expr*), tanh (*expr*), asinh (*expr*), acosh (*expr*), atanh (*expr*)
    *expr-valued:* The hyperbolic sine, cosine, tangent of *expr* and the corresponding
    inverse functions.
- *expr1* ^ *expr2*
    *expr-valued:* *expr1* raised to the *expr2* power, for the special cases where
    either *expr1* or *expr2* is a constant. For *expr2* positive integer, the operator
    is decomposed into quadratic constraints if the solver supports them,
    otherwise passed to the solver natively or approximated by a piecewise-linear function.

General solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For linear-quadratic MP-based solvers,
most of these nonlinear functions are handled by piecewise-linear approximation,
except products with binary variables.
The appoximation is constructed by the MP interface, using options
*cvt:plapprox:reltol* and *cvt:plapprox:domain*,
and is then processed as described in
:ref:`piecewise_linear_modeling`.

Handling in Gurobi
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For Gurobi, the following univariate nonlinear functions are instead handled natively:
**exp**, **log**, **^**, **sin**, **cos**, **tan**.
After suitable transformations, the MP interface sends Gurobi the expressions that use
these functions.

Gurobi 11 defaults to piecewise-linear approximation of these functions
as part of its preprocessing. However, Gurobi :ref:`option <solver-options>` *alg:global*
(*pre:funcnonlinear*) sets the default to the new MINLP capability ---
`global nonlinear solving <https://www.gurobi.com/>`_ via spatial branching::

  ampl: option gurobi_options 'global=1'; solve;

For individual constraints and objectives, the choice of global solving vs
piecewise-linear approximation can be performed via the `.global` suffix::

  pre:funcnonlinear (funcnonlinear, global)
      Controls how general functions with their constraint's or objective's
      suffix .funcnonlinear or, if not available, .global unset (or set to 0)
      are treated (ATTENTION: different meaning than Gurobi FuncNonLinear
      parameter and attribute):

      -1 - Piecewise-linear approximation
      0  - Automatic (default)
      1  - Treated as nonlinear functions

      Suffix values mean the same.

See our `Christmas tree decorations Streamlit example <https://ampl.com/streamlit>`_
using global optimization.

Piecewise-linear aproximation can be influenced by setting
the following options in the AMPL ``gurobi_options`` string::

  pre:funcpieces
      Sets the strategy for constructing a piecewise-linear approximation of a
      function:

      0   - Automatic choice (default)
      >=2 - Sets the number of pieces, of equal width
      1   - Uses a fixed width for each piece, as specified by the
            funcpiecelength option
      -1  - Bounds the absolute error of the approximation, as specified
            by the funcpieceerror option
      -2  - Bounds the relative error of the approximation, as specified
            by the funcpieceerror option

  pre:funcpiecelength
      When funcpieces = 1, specifies the length of each piece of the
      approximation.

  pre:funcpieceerror
      When funcpieces = -1 or -2, specifies the maximum allowed
      error (absolute for -1, relative for -2) in the approximation.

  pre:funcpieceratio
      Controls whether the piecewise-linear approximation is an underestimate
      of the function, an overestimate, or somewhere in between. A value of
      0.0 will always underestimate, while a value of 1.0 will always
      overestimate; a value in between will interpolate between the
      underestimate and the overestimate. A special value of -1 chooses
      points that are on the original function.

These options can also be overridden for a particular objective or constraint,
by setting suffixes of the same names. For example, after defining the objective
shown below, setting ``suffix funcpieces IN; let Chichinadze.funcpieces := 12;``
specifies 12 pieces for approximating the sin, cos, and exp functions in that objective.

.. code-block:: ampl

    minimize Chichinadze:
       x[1]^2 - 12*x[1] + 11 + 10*cos(pi*x[1]/2) +
          8*sin(pi*5*x[1]) - exp(-(x[2]-.5)^2/2)/sqrt(5);

