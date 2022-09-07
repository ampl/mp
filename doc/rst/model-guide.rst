.. _modeling-guide:

Modeling Guide 
=========================
for MP-based AMPL Solvers
-------------------------

AMPL's newly extended C++ solver interface library, MP, is publicly available in the `ampl/mp <https://github.com/ampl/mp>`_ repository. Solver interfaces built with MP are able to handle a significantly expanded range of model expressions. Currently available MP-based solvers include:

- `x-gurobi <https://github.com/ampl/mp/tree/master/solvers/gurobi>`_, an enhanced interface to the `Gurobi <https://ampl.com/products/solvers/solvers-we-sell/gurobi/>`_ solver

- `copt <https://github.com/ampl/mp/tree/master/solvers/copt>`_, an interface to `Cardinal Optimizer <https://ampl.com/products/solvers/solvers-we-sell/copt/>`_

- `highs <https://github.com/ampl/mp/tree/master/solvers/highsdirect>`_, an interface to the open-source `HiGHS <https://highs.dev/>`_ solver

Binaries for these solvers can be downloaded, in distribution bundles and individually, through the `AMPL Portal <https://portal.ampl.com>`_.


Overview
--------

The expanded MP solver interface library offers new support for the following categories of operators and expressions:

- Conditional operators: ``if-then-else``; ``==>``, ``<==``, ``<==>``
- Logical operators: ``or``, ``and``, ``not``; ``exists``, ``forall``
- Piecewise linear functions: ``abs``; ``min``, ``max``; ``<<breakpoints; slopes>>``
- Counting operators: ``count``; ``atmost``, ``atleast``, ``exactly``; ``numberof``
- Comparison operators: ``>``, ``<``, ``!=``; ``alldiff``
- Complementarity operator: ``complements``
- Nonlinear operators and functions: ``*``, ``/``, ``^``; ``exp``, ``log``; ``sin``, ``cos``, ``tan``
- Set membership operator: ``in``

Details and examples are given in the *Expressions supported* section below. See also the individual solvers' documentation for details of solver-specific features:

- Choice between linearization in the interface and native solver support for some operations
- Handling of AMPL suffixes on constraints that are transformed by the interface

The slides from our presentation on `Advances in Model-Based Optimization <https://ampl.com/MEETINGS/TALKS/2022_07_Bethlehem_Fourer.pdf>`_ provide overview of the MP interface library in the context of AMPL applications, including comments on implementation and efficiency issues. 


Expressions supported
---------------------

The MP solver interface library works with existing AMPL syntaxes, but allows them to be used in more general ways, or with a greater variety of solvers.

In many cases, an extension results from allowing variables to appear in more general contexts, such as with conditional, logical, or counting operators. Other extensions are enabled by providing more powerful transformations, particularly to linear or quadratic equivalents, and by providing support for extensions that are native to some solvers. A few extensions are already handled in the AMPL language translator, and are included here for completeness.

In the syntax summaries below, there are two main kinds of entities, representing *numerical expressions* and *constraints:*

- **expr** 
     represents any expression that evaluates to a number. Unless otherwise indicated, it may contain variables. It may be built from familiar arithmetic operators, but also from other operators or functions that return numerical values.

- **constr** 
     represents a constraint of the model, which may evaluate to true or false depending the values of variables that it contains. It may be built from the familiar relational operators ``>=``, ``<=``, and ``=``, but also from other operators such as ``or`` and ``alldiff`` that create constraints.

The return value of an operator or function is also either an *expr* or *constr*, as indicated. Thus it is possible to build up complex combinations of numerical and logical operators; for example,

.. code-block:: ampl

        (x<=0 or y!=2)  ==>
                (x<=-5 or
                        (max((x+1)*(x+2)*(y+3), y)<=3 and exp((x+18)*y)<=12));

AMPL represents these as expression trees, which are sent to MP-based solver interfaces to be processed as particular solvers require.


Conditional operators
***********************************

- if *constr* then *expr1* [else *expr2*]
    *Returns expr:* When *constr* is true, takes the value of *expr1*.  
    When *constr* is false, takes the value of *expr2*, or 0 if the `else` phrase is omitted.

In the special case where there are no variables in the *constr*, the value of this expression can be determined as either *expr1* or *expr2* (or 0) before the problem is sent to the solver. But in general, the expression's value depends upon how the solver sets the variables in the *constr*, and so AMPL must send the entire expression to the solver interface for processing.

.. code-block:: ampl

       minimize TotalCost:
          sum {j in JOBS, k in MACHINES}
             if MachineForJob[j] = k then cost[j,k];

.. code-block:: ampl

       subject to Balance {p in PROD, t in TIME}:
          Make[p,t] + (if t = 0 then inv0[p] else Inv[p,t-1])
             = Sell[p,t] + Inv[p,t];

- *constr1* ==> *constr2*
- *constr2* <== *constr1*
    *Returns constr:* Satistifed if *constr1* is true and *constr2* is true, 
    or if *constr1* is false. 
- *constr1* ==> *constr2* else *constr3*
    *Returns constr:* Satistifed if *constr1* is true and *constr2* is true, 
    or if *constr1* is false and *constr3* is true. 
- *constr1* <==> *constr2*
    *Returns constr:* Satisfied if *constr1* and *constr2* are both true or both false.


Logical operators
***********************************

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


General combinatorial expressions
*********************************

SOS constraints and non-contiguous variable domains
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SOS1 is mainly relevant for models that restrict some variables to take a
value from an arbitrary list of values. A simple example:

.. code-block:: ampl

    var Buy {f in FOODS} in {0,10,30,45,55};

An appropriate SOS1 representation will be
automatically generated from this declaration.

SOS2 are one of the two ways to linearize
:ref:`piecewise-linear expressions <piecewize-linear-expr>` by AMPL.

It is possible to specify SOS1 or SOS2 variables and corresponding "reference rows"
explicitly using AMPL suffixes .sos(no) and .(sos)ref,
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


.. _piecewize-linear-expr:

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



Nonlinear expressions
*********************


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

The piecewise-linear approximation is controlled by :ref:`Gurobi-FuncPieces`.


Suffix conversions
------------------

MP converts suffixes between the original and transformed model
('value presolve'), in particular *irreducible independent subsystem* (IIS)
results and Gurobi `FuncPieces` and related attributes.


IIS reporting
*************

As an example, for the following model:

.. code-block:: ampl

    var x;
    var y;
    var z;

    subj to Con1:
       x+y >= 1;

    subj to Con2:
       y + log(z + exp(x+3)) <= 1.83;

    subj to Con3:
       z + log(y + 3.8*exp(x+3)) >= -14.265;

all constraints are reported as IIS members:

.. code-block:: ampl

    ampl: option gurobi_options 'iisfind=1';
    ampl: solve;
    ....
    ampl: display _con.iis;
    _con.iis [*] :=
    1  mem
    2  mem
    3  mem
    ;


.. _Gurobi-FuncPieces:

Gurobi `FuncPieces` and related parameters
******************************************

Gurobi functional constraint attributes `FuncPieces`, `FuncPieceLength`,
`FuncPieceError`, and `FuncPieceRatio` determine the piecewise-linear
approximation applied. The MP Gurobi driver defines the corresponding
options relating to the whole model, but also suffixes for constraints,
which are converted to Gurobi representation. Example: for the above
IIS model, setting the `.funcpieces` suffix as follows:

.. code-block:: ampl

    suffix funcpieces IN;

    let Con1.funcpieces := 12;
    let Con2.funcpieces := 23;
    let Con3.funcpieces := 38;

results in the following Gurobi model (LP format, excerpt):

.. code-block:: ampl

    ...
    General Constraints
     GC0: ( FuncPieces=38 ) C4 = EXP ( C3 )
     GC1: ( FuncPieces=23 ) C6 = LOG ( C5 )
     GC2: ( FuncPieces=38 ) C8 = LOG ( C7 )
    End



Conversion graph export
-----------------------

The conversion graph can be exported using the `writegraph` option,
currently in JSON Lines format.


Efficient modeling
------------------

For general modeling advice, refer to sources such as
Guidelines for Numerical Issues
and modeling webinars on the `Gurobi website <http://www.gurobi.com>`_,
Practical Considerations for Integer Programming in the
`AMPL Book <https://ampl.com/resources/the-ampl-book/>`_, and
the MOSEK Modeling Cookbook at `www.mosek.com <https://www.mosek.com/>`_.


Reduce non-linearity
********************

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
************

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
