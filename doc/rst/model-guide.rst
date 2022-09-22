.. _modeling-guide:

Modeling Guide 
=========================
for MP-based AMPL Solvers
-------------------------

AMPL's newly extended C++ solver interface library, MP, is publicly available in the `ampl/mp <https://github.com/ampl/mp>`_ repository. Solver interfaces built with MP are able to handle a significantly expanded range of model expressions. Currently available MP-based solvers include:

- `x-gurobi <https://github.com/ampl/mp/tree/master/solvers/gurobi>`_, an enhanced interface to the `Gurobi <https://ampl.com/products/solvers/solvers-we-sell/gurobi/>`_ solver

- `copt <https://github.com/ampl/mp/tree/master/solvers/copt>`_, an interface to `Cardinal Optimizer <https://ampl.com/products/solvers/solvers-we-sell/copt/>`_

- `highs <https://github.com/ampl/mp/tree/master/solvers/highsdirect>`_, an interface to the open-source `HiGHS <https://highs.dev/>`_ solver

Binaries for these solvers can be downloaded, in distribution bundles and individually, through the `AMPL Portal <https://portal.ampl.com>`_. More solvers will be added.


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

The return value of an operator or function is also one of the above, as indicated by *expr-valued* or *constr-valued* at the beginning of each syntax summary. Thus it is possible to build up complex combinations of operators and functions of various kinds; for example,

.. code-block:: ampl

        (x<=0 or y!=2)  ==>
                (x<=-5 or
                        (max((x+1)*(x+2)*(y+3), y)<=3 and exp((x+18)*y)<=12));

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
        max (cost["BRO"],cost["CAU"],cost["BRU")
        max ({f in FOOD} cost[f], 10.0)
        max ({n in NUTR} (lim_nutr[n], {f in FOOD} amt[n,f])) 

As seen in the case of ``max`` above, certain operators can be used with either the ``{indexing} expr`` or the ``(expr-list)`` form.

Due to the generality of the operators recognized by the MP interface, it is possible to express constraints that do not define a closed feasible region. For example,

.. code-block:: ampl

        x > 5
        not (x >= 5 and y >= 15)
        x = 0 ==> z = 0 else z = 1

An optimum is not guaranteed to exist over a non-closed region. Thus where necessary, the MP interface constructs an approximate closed region by use of a small tolerance. For example, if x is minimized subject to x > 5, then any x value greater than 5 is not minimal, and any x value less than or equal to 5 is not feasible. Thus, to insure that a minimum is well defined, the constraint must be changed to x >= 5 + eps for some small constant eps. Each solver has its own default value of the eps constant, which can be adjusted through an option setting. 


Conditional operators
***********************************

- if *constr* then *expr1* [else *expr2*]
    *expr-valued:* When *constr* is true, takes the value of *expr1*.  
    When *constr* is false, takes the value of *expr2*, or 0 if the `else` phrase is omitted.

In the special case where there are no variables in the *constr*, the value of this expression can be determined as either *expr1* or *expr2* (or 0) before the problem is sent to the solver. But in general, the value of expression depends upon how the solver assigns values to the variables in the *constr*, and so AMPL must send the entire expression to the solver interface for processing.

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

The conditional expression *constr1* ==> *constr2* can be thought of as saying that *constr1* implies *constr2*, or equivalently that if *constr1* then *constr2*. In the special case where *constr1* is of the form *binary-var* = 0 or *binary-var* = 1, these are "indicator" constraints that can be handled natively by some solvers. Otherwise, they are transformed to simpler constraints that use relational operators. The other cases are treated similarly.

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
    
Expressions using these operators are transformed to use Gurobi's native AND and OR "general constraints" when possible. In other cases, they are transformed to simpler constraints that use relational operators.

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
    *expr-valued:* Computes a piecewise-linear function of a single variable; see `Chapter 17. Piecewise-Linear Programs <https://ampl.com/BOOK/CHAPTERS/20-piecewise.pdf>`_ in the `AMPL book <https://ampl.com/resources/the-ampl-book/>`_ for a complete description of the forms that AMPL recognizes.
    
This piecewise-linear expression is defined by lists of ``n`` *breakpoints* and ``n+1`` *slopes*. The *var* must be a reference to a single variable.

When AMPL's option ``pl_linearize`` is at its default value of 1, AMPL linearizes these piecewise-linear expressions, and sends the linearized versions to the solver. The linearization is continuous where possible, in certain convex and concave cases (where the slopes are increasing and decreasing, respectively); but in general, the linearization includes both continuous and binary variables.

When ``pl_linearize`` is set to 0, piecewise-linear expressions are represented to the solver in the form of expression trees. The MP-based interface transforms them to use Gurobi's native methods for piecewise-linear functions, and linearizes them for other solvers.

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
    *expr-valued:* The number of items in the *expr-list* have the same value as *expr*.

This operator can provide an easier-to-read alternative for a special case of count. Compare for example the ``CapacityOfMachine`` constraint below to the one given previously using ``atmost``.

.. code-block:: ampl

    subj to CapacityOfMachine {k in MACHINES}:
        numberof k in ({j in JOBS} MachineForJob[j]) <= cap[k];

.. code-block:: ampl

    subj to MinInGrpDefn {j in 1..numberGrps}:  
       MinInGrp <= numberof j in ({i in PEOPLE} Assign[i]);


Comparison operators
***********************************

- expr1 > expr2
    *constr-valued:* Satisfied when *expr1* is strictly greater than *expr2*.
- expr1 < expr2
    *constr-valued:* Satisfied when *expr1* is strictly less than *expr2*.
- expr1 != expr2
    *constr-valued:* Satisfied when *expr1* does not equal *expr2*.

Where possible, the MP interface transforms these operations to ones involving ``>=`` and ``<=``, so that optimization solvers can handle them. For example, this can be done when *expr1* and *expr2* are integer-valued, or when an expression like ``if Flow[i,j] > 0 then fixed[i,j]`` expresses a fixed cost in an objective to be minimized. Where this is not possible, a small tolerance is introduced as disucssed in the section above on **Expressions supported**.

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

This operator provides a much more concise alternative to specifying ``!=`` between all pairs in a specified collection of expressions. Currently none of the MP-based solvers support this operator natively, so the interface transforms it to a representation in terms of simpler constraints that use relational operators.

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
- **expr complements constr**
- **constr complements expr**
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
