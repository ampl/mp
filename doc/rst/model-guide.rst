.. _modeling_guide:

Modeling guide
==============

A guide to using the beta test release of
`x-gurobi <https://github.com/ampl/mp/tree/master/solvers/gurobidirect>`_,
the enhanced
`AMPL-Gurobi <https://ampl.com/products/solvers/solvers-we-sell/gurobi/>`_
interface, and
`copt <https://github.com/ampl/mp/tree/master/solvers/copt>`_, an interface
to `Cardinal Optimizer <https://www.shanshu.ai/copt>`_.
Both
can be compiled from source or downloaded in the AMPL distribution bundle.


Summary
-------

- Full support of logical expressions and constraints, as described in the
  AMPL page on `Logic and Constraint Programming Extensions
  <https://ampl.com/resources/logic-and-constraint-programming-extensions/>`_.
  
- Algebraic expressions beyond linear and quadratic, including
  complementarity constraints.

- Choice between conversions in the driver vs. native solver support.


Expressions supported
---------------------

- Arbitrary trees of logical, relational, and general non-linear expressions
  including higher-degree polynomials:

  .. code-block:: ampl

        (x<=0 or y!=2)  ==>
                (x<=-5 or
                        (max((x+1)*(x+2)*(y+3), y)<=3 and
                                exp((x+18)*y)<=12));

- Piecewise-linear expressions are natively supported by Gurobi.
  To pass them that way, set the following AMPL option:

  .. code-block:: ampl

        option pl_linearize 0;

- Complementarity constraints:

  .. code-block:: ampl

        subject to Pri_Compl {i in PROD}:
            max(500.0, Price[i]) >= 0 complements
                sum {j in ACT} io[i,j] * Level[j] >= demand[i];

        subject to Lev_Compl {j in ACT}:
            level_min[j] <= Level[j] <= level_max[j] complements
                cost[j] - sum {i in PROD} Price[i] * io[i,j];

- Constraint programming high-level constraints and expressions, for example:

  .. code-block:: ampl

        ## if-then
        minimize TotalCost:
            sum {j in JOBS, k in MACHINES}
                if MachineForJob[j] = k then cost[j,k];

        ## implied forall
        subj to HostNever {j in BOATS}:
            isH[j] = 1 ==> forall {t in TIMES} H[j,t] = j;

        ## alldiff global constraint
        subj to OneJobPerMachine:
            alldiff {j in JOBS} MachineForJob[j];

        ## implied alldiff
        subj to VisitOnce {j in BOATS}:
            isH[j] = 0 ==> alldiff {t in TIMES} H[j,t];

        ## numberof operator
        subj to CapacityOfMachine {k in MACHINES}:
            numberof k in ({j in JOBS} MachineForJob[j]) <= cap[k];

        ## count operator
        subj to CapacityOfMachine {k in MACHINES}:
            count {j in JOBS} (MachineForJob[j] = k) <= cap[k];

        ## implied atmost
        subj to VisitHosts {i in BOATS}:
            isH[i] = 0 ==> atmost 0 {j in BOATS, t in TIMES} (H[j,t] = i);


- QP expressions are multiplied out. Example:

  .. code-block:: ampl

        -5 * (abs(x[1])-0.7)^2 + x[2]

  is converted as follows:

  .. code-block:: ampl

      -5*t^2 + 7*t - 2.45 + x[2]

      s.t. t = abs(x[1]);

- Gurobi “general constraints” ``and``, ``or``, ``max``, ``min``, ``abs``,
  as well as indicator (``==>``), are passed to Gurobi natively.
  This behavior can be changed with solver options **acc:abs** etc.
  to use big-*M* constraints instead (when the variables have
  finite bounds).

- Nonlinear “generals” are passed to Gurobi:

  ``exp``, ``log``, ``sin``, ``cos``, ``tan``, ``pow``, ``pl``, ``SOS1/2``.

General modeling hints
----------------------

For general modeling hints, refer to Guidelines for Numerical Issues
and modeling webinars on the `Gurobi website <http://www.gurobi.com>`_;
Practical Considerations for Integer Programming in the
`AMPL Book <https://ampl.com/resources/the-ampl-book/>`_, and
the MOSEK Modeling Cookbook at `www.mosek.com <https://www.mosek.com/>`_.

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
