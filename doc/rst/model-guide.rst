.. _modeling_guide:

x-Gurobi modeling guide
=======================

A guide to using the beta test release of
`x-gurobi <https://github.com/ampl/mp/tree/master/solvers/gurobidirect>`_,
the enhanced AMPL-Gurobi interface.
`x-gurobi <https://github.com/ampl/mp/tree/master/solvers/gurobidirect>`_
can be compiled from source or downloaded in the AMPL distribution bundle.


Summary
-------

- Full support of logical expressions and constraints, as described in the
  AMPL page on `Logic and Constraint Programming Extensions
  <https://ampl.com/resources/logic-and-constraint-programming-extensions/>`_.
  
- Algebraic expressions beyond linear and quadratic.

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


- QP expressions are multiplied out:

  .. code-block:: ampl

        -5 * (abs(x[1])-0.7)^2 + x[2]

- Gurobi “general constraints” ``and``, ``or``, ``max``, ``min``, ``abs``,
  as well as indicator (``==>``), are passed to Gurobi natively.
  This behavior can be changed with solver options **acc:abs** etc.
  to use big-*M* constraints instead (when the variables have
  finite bounds).

- Nonlinear “generals” are passed to Gurobi:

  ``exp``, ``log``, ``sin``, ``cos``, ``tan``, ``pow``, ``pl``, ``SOS1/2``.

General modeling hints
----------------------

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
