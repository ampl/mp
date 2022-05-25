.. _modeling-guide:

Modeling guide
==============

A guide to using the beta test release of
`x-gurobi <https://github.com/ampl/mp/tree/master/solvers/gurobi>`_,
the enhanced
`AMPL-Gurobi <https://ampl.com/products/solvers/solvers-we-sell/gurobi/>`_
interface,
`copt <https://github.com/ampl/mp/tree/master/solvers/copt>`_, an interface
to `Cardinal Optimizer <https://www.shanshu.ai/copt>`_, and
`highs <https://github.com/ampl/mp/tree/master/solvers/highsdirect>`_, an interface
to `HiGHS <https://highs.dev/>`_.
They can be downloaded in the `AMPL distribution bundle <https://portal.ampl.com>`_
or compiled from source.


Summary
-------

- Full support of logical expressions and constraints, as described in the
  AMPL page on `Logic and Constraint Programming Extensions
  <https://ampl.com/resources/logic-and-constraint-programming-extensions/>`_.
  
- Algebraic expressions beyond linear and quadratic, including
  complementarity constraints.

- Choice between conversions in the driver vs. native solver support.

Example models can be found in the
`test suite <https://github.com/ampl/mp/tree/develop/test/end2end/cases>`_ and
in the
`ampl.github.io model collection <https://github.com/ampl/ampl.github.io/tree/master/models>`_.
See also an `overview talk <https://ampl.com/MEETINGS/TALKS/2022_04_Houston_Tutorial.pdf>`_.


Expressions supported
---------------------

MP supports arbitrary trees of logical, relational, and general non-linear expressions
including higher-degree polynomials:

.. code-block:: ampl

        (x<=0 or y!=2)  ==>
                (x<=-5 or
                        (max((x+1)*(x+2)*(y+3), y)<=3 and exp((x+18)*y)<=12));

Below are details on the various kinds of expressions and how they are presented
to the solvers.


Conditional expressions
~~~~~~~~~~~~~~~~~~~~~~~

AMPL already provides an ``if-then-else`` operator that returns a value
that can be used in expressions:

- if *logical-expr* then *object-expr1*

- if *logical-expr* then *object-expr1* else *object-expr2*
    Takes the value of object-expr1 when the *logical-expr* is true, and the value
    of *object-expr2* (or 0 if omitted) when the *logical-expr* is false.
    Both *object-expr*'s must be compatible, i.e., both numbers, strings, or sets.


When this operator appears in a constraint, the *logical-expr*
can contain variables, in which case AMPL handles the constraint like
other nonlinear constraints, passing an expression tree to the solver.
In particular, the *logical-expr* may be any valid *constraint-expr*.

.. code-block:: ampl

        ## if-then
        minimize TotalCost:
            sum {j in JOBS, k in MACHINES}
                if MachineForJob[j] = k then cost[j,k];


Logical relations
~~~~~~~~~~~~~~~~~

AMPL also has a similar if-then form of indexing expression,
which is used in the context of constraints as follows:

- subject to name {if *logical-expr*}: *constraint-expr*;
    Enforces the constraint specified by the *constraint-expr*
    if and only if the *logical-expr* is true.

Thus for example in section 8.4 of the
`AMPL book <https://ampl.com/resources/the-ampl-book/>`_ we have:

.. code-block:: ampl

     subject to Time {if avail > 0}:
         sum {p in PROD} (1/rate[p]) * Make[p] <= avail;

It is arguably more natural, however, to make the ``if`` condition part of the
constraint expression. Since the ``if-then`` and ``if-then-else`` constructs
are already heavily used in AMPL (for expressions and for script statements),
we have introduced several operators for describing implications in constraints.
For example:

.. code-block:: ampl

    subject to Time:
        avail > 0 ==> sum {p in PROD} (1/rate[p]) * Make[p] <= avail;

General forms of AMPL’s logical relations are as follows:

- *logical-expr* ==> *constraint-expr1*
    Satisfied if the *logical-expr* is true and *constraint-expr1* is satisfied,
    or if the *logical-expr* is false.
- *logical-expr* ==> *constraint-expr1* else *constraint-expr2*
    Satisfied if the *logical-expr* is true and *constraint-expr1* is satisfied,
    or if the *logical-expr* is false and *constraint-expr2* is satisfied.
- *logical-expr* <==> *constraint-expr*
    Satisfied if the *logical-expr* is true and *constraint-expr* is satisfied,
    or if the *logical-expr* is false and *constraint-expr* is not satisfied.

Additionally ``<==`` has the same meaning as ``==>`` except with the roles of
*constraint-expr1* and *constraint-expr2* reversed.

By allowing variables on both sides of the implication operators,
these forms considerably expand the variety of conditional constraints
that AMPL can conveniently express. For example:

.. code-block:: ampl

    subject to Multi_Min_Ship {i in ORIG, j in DEST}:
        sum {p in PROD} Trans[i,j,p] > 0 ==>
            minload <= sum {p in PROD} Trans[i,j,p] <= limit[i,j];

Again, the *logical-expr* can be any *constraint-expr*.
Conditional operators can be nested and combined with other operators.

AMPL conditional operators are either linearized using big-*M* constraints, or passed
to the solver natively as indicator constraints
(if supported; e.g., Gurobi options *acc:ind_le*, *acc:ind_eq*).


Logical expressions
~~~~~~~~~~~~~~~~~~~

Basic AMPL constraints consist of numerical-valued expressions
connected by ``<=``, ``>=`` or ``=``. These constraint expressions
are now allowed to be
connected by AMPL’s unary and binary logical operators,

- *constraint-expr1* or *constraint-expr2*
    Satisfied iff at least one of the operands is satisfied.
- *constraint-expr1* and *constraint-expr2*
    Satisfied iff both of the operands are satisfied.
- not *constraint-expr*
    Satisfied iff the operand is not satisfied.

and AMPL’s iterated forms of the binary logical operators:

- exists {indexing} *constraint-expr*
    Satisfied iff the operand is satisfied for at least one
    member of the indexing set (the iterated form of ``or``).
- forall {indexing} *constraint-expr*
    Satisfied iff the operand is satisfied for all members of
    the indexing set (the iterated form of ``and``).
- forall ( {indexing} *constraint-expr1*, {indexing} *constraint-expr2*, ...)
    Example of compound indexing. Each {indexing} may be any AMPL
    indexing-expression, or may be omitted to specify a single
    item in the list.

.. Meaning of the below?
  Constraint expressions can also be grouped by parentheses:
  ( constraint-expr )
  Satisfied iff the constraint-expr is satisfied.

So an AMPL constraint can be any logical combination of equalities,
inequalities and other boolean expressions:

.. code-block:: ampl

        subj to HostNever {j in BOATS}:
            isH[j] = 1 ==> forall {t in TIMES} H[j,t] = j;

Using the ``not`` operator it is possible to specify a feasible region
that isn’t closed, so that optimization problems using continuous
variables may be meaningless. This is illustrated by a very simple problem:

.. code-block:: ampl

    var x;
    minimize Obj: x;
    subject to OpenCons: not (x <= 2);

The objective has an infimum of 2, but no minimum that satisfies the
constraint. The same problem arises if one uses a strict inequality ``<``
or ``>``, specifically the expresion ``x > 2`` in this case.
For MIP solvers, MP redefines strict inequalities using a tolerance
(option *cvt:mip:eps*).
Most CP solvers, operating only on discrete variables,
freely allow expressions that have these forms.


AMPL logical expressions are either linearized using boolean arithmetic, or passed
to the solver natively
(if supported; e.g., Gurobi options *acc:and*, *acc:or*).


SOS variable domains
~~~~~~~~~~~~~~~~~~~~

SOS1 is mainly relevant for models that restrict some variables to take a
value from an arbitrary list of values. A simple example:

.. code-block:: ampl

    var Buy {f in FOODS} in {0,10,30,45,55};

An appropriate SOS1 representation will be
automatically generated from this declaration.

It is possible to specify SOS1 or SOS2 variables and corresponding "reference rows"
explicitly using AMPL suffixes .sosno and .ref,
as described in the solver documentation.
However this requires some study to understand whether SOS1/2 is appropriate
and how to apply it, and we don't recommend going to that trouble unless you
are having serious problems getting the solver to return a solution.


Min, max, abs
~~~~~~~~~~~~~

Non-smooth functions ``min`` and ``max`` can have either a fixed argument list,
or be iterated:

.. code-block:: ampl

    abs(x)
    min(x, y, max(z, 2))
    max {i in ORIG} supply[i]

Functions ``min``, ``max``, ``abs`` can be linearized with big-*M* constraints
or passed to the solver natively
(if supported; e.g., Gurobi options *acc:min*, *acc:max*, *acc:abs*).


Piecewise-linear expressions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A piecewise-linear expression is defined by a list of ``n`` *breakpoints*
and ``n+1`` *slopes*, together with an argument variable:

.. code-block:: ampl

    <<limit1[i,j], limit2[i,j];
      rate1[i,j], rate2[i,j], rate3[i,j]>> Trans[i,j]

In this example, ``n=2`` and the argument is the variable ``Trans[i,j]``.
An AMPL PL expression
assumes that the corresponding function passes through origin (0, 0).
See the `AMPL book <https://ampl.com/resources/the-ampl-book/>`_
for more information.

Solvers natively supporting piecewise-linear expressions,
for example Gurobi, perform best when receive them that way
(vs linearization by AMPL, which is currently the default).
To do so, switch off the corresponding AMPL option:

.. code-block:: ampl

        option pl_linearize 0;


Complementarity constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~

AMPL accepts two kinds of complementarity constraints.
The first kind, inequality vs inequality, enforces both inequalities
and makes sure at least one of them is tight:

.. code-block:: ampl

        subject to Pri_Compl {i in PROD}:
            max(500.0, Price[i]) >= 0 complements
                sum {j in ACT} io[i,j] * Level[j] >= demand[i];

The second kind, range constraint vs expression,
enforces one of the following 3 cases:

1. range constraint at lower bound  and  expression >= 0;
2. range constraint valid and expression == 0;
3. range constraint at upper bound and expression <= 0, for example:

.. code-block:: ampl

        subject to Lev_Compl {j in ACT}:
            level_min[j] <= Level[j] <= level_max[j] complements
                cost[j] - sum {i in PROD} Price[i] * io[i,j];

See the `AMPL book <https://ampl.com/resources/the-ampl-book/>`_
for more information.

Quadratic expressions are allowed. For MIP solvers, complementarity
conditions are represented by logical constraints.


Counting operators
~~~~~~~~~~~~~~~~~~

AMPL’s ``count`` operator returns the number of times that
a certain constraint is satisfied:

- count {indexing} *constraint-expr*
    The number of members of the indexing set such that the
    *constraint-expr* is satisfied.

The *constraint-expr* can be any valid AMPL constraint.
The AMPL translator will instantiate it for each member of
the indexing set, and will communicate all of the instantiated
constraints to the solver interface.

Additional iterated logical operators are provided to simplify
the descriptions of constraints in some common special cases:

- atmost k {indexing} *constraint-expr*
    Satisfied iff the *constraint-expr* holds for at most ``k`` members of the indexing set.
- atleast k {indexing} *constraint-expr*
    Satisfied iff the *constraint-expr* holds for at least ``k`` members of the indexing set.
- exactly k {indexing} *constraint-expr*
    Satisfied iff the *constraint-expr* holds for exactly ``k`` members of the indexing set.

``k`` can be any constant arithmetic expression that evaluates to a nonnegative integer value.

Another particularly important special case occurs when counting the number of set members
at which a given expression takes a particular value.
The general form is:

- numberof k in ({indexing} *object-expr*)
    The number of members of the indexing set such that the *object-expr* is equal to ``k``.


.. code-block:: ampl

        ## numberof operator
        subj to CapacityOfMachine {k in MACHINES}:
            numberof k in ({j in JOBS} MachineForJob[j]) <= cap[k];

        ## implied atmost
        subj to VisitHosts {i in BOATS}:
            isH[i] = 0 ==> atmost 0 {j in BOATS, t in TIMES} (H[j,t] = i);


Pairwise operator
~~~~~~~~~~~~~~~~~

Various assignment and related combinatorial problems require that
a collection of entities be pairwise different or disjoint. Operator ``alldiff``
makes these conditions easier to state and helps to make the resulting problems
easier to solve.

In general, this operator can be applied to any collection of expressions
involving variables:

- alldiff {indexing} *var-expr*
- alldiff ( {indexing} *var-expr1*, {indexing} *var-expr2*, ... )
    Satisfied iff all of the specified variables take different values. Each
    {indexing} may be any AMPL indexing-expression, or may be omitted to
    specify a single item in the list.

.. code-block:: ampl

        ## implied alldiff
        subj to VisitOnce {j in BOATS}:
            isH[j] = 0 ==> alldiff {t in TIMES} H[j,t];


QP and polynomials
~~~~~~~~~~~~~~~~~~

QP expressions are multiplied out. For example, the following expression:

.. code-block:: ampl

    -5 * (abs(x[1])-0.7)^2 + x[2]

is converted as follows:

.. code-block:: ampl

    -5*t^2 + 7*t - 2.45 + x[2]

with an auxiliary variable ``t = abs(x[1])``.

Higher-order algebraic expressions are broken down to quadratics
via auxiliary variables:

.. code-block:: ampl

    maximize Sum:
        -5 * (x[1]-0.7)^2 + x[2]^7;


Nonlinear functions
~~~~~~~~~~~~~~~~~~~

Gurobi 9 introduced non-linear functional constraints which are internally
handled by piecewise-linear approximation. The following are the corresponding
AMPL functions:

``exp``, ``log``, ``sin``, ``cos``, ``tan``, ``pow``.



Efficient modeling
------------------

For general modeling advice, refer to Guidelines for Numerical Issues
and modeling webinars on the `Gurobi website <http://www.gurobi.com>`_,
Practical Considerations for Integer Programming in the
`AMPL Book <https://ampl.com/resources/the-ampl-book/>`_, and
the MOSEK Modeling Cookbook at `www.mosek.com <https://www.mosek.com/>`_.


Reduce non-linearity
~~~~~~~~~~~~~~~~~~~~

In the following example:

.. code-block:: ampl

    var Flow {PRODUCTS,ARCS} >= 0;

    minimize TotalCost:
        sum {(i,j) in ARCS}
            if exists {p in PRODUCTS} Flow[p,i,j] > 0 then fix_cost[i,j];

it is possible to reduce the number of resulting indicator constraints
via the following simplification:

.. code-block:: ampl

    minimize TotalCost:
        sum {(i,j) in ARCS}
            if sum {p in PRODUCTS} Flow[p,i,j] > 0 then fix_cost[i,j];

Such a simplification might be performed automatically in a future version
of the library.


Tight bounds
~~~~~~~~~~~~

For logical expressions, it proves best to supply tight bounds on
all participating variables.
For any intermediate expressions which are known to have tighter bounds
than those which can be deduced automatically, it is advisable
to extract them into extra variables with the tight bounds.
For example, given a disjunction

.. code-block:: ampl

        subj to: log(x+2)<=y^2  or  x-y>=z;

and knowing that ``-15 <= x-y-z <= 30``, reformulate:

.. code-block:: ampl

        var t >=-15, <=30;
        subj to: t == x-y-z;
        subj to: log(x+2)<=y^2  or  t>=0;

In many cases, integer variables are more meaningful and efficient
in logical constraints
than continuous variables, for example in disequalities.
