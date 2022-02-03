x-Gurobi Modeling Guide
=======================

A guide to using the beta test release of
`x-gurobi <https://github.com/ampl/mp/tree/master/solvers/gurobidirect>`_,
the enhanced AMPL-Gurobi interface.


Summary
-------

- Full support of logical expressions and constraints, as described in the
  AMPL page on `Logic and Constraint Programming Extensions
  <https://ampl.com/resources/logic-and-constraint-programming-extensions/>`_.
  
- Algebraic expressions beyond linear and quadratic.

- Choice between conversions in the driver vs. native solver support.


Expressions supported
---------------------

- Logical, relational, and general including higher-degree polynomials:

  .. code-block:: ampl

        (x<=0 or y!=2)  ==>
                (x<=-5 or
                        (max((x+1)*(x+2)*(y+3), y)<=3 and b==0));

- Constraint programming high-level constraints and expressions, for example:

  .. code-block:: ampl

        minimize TotalCost:
            sum {j in JOBS, k in MACHINES}
                if MachineForJob[j] = k then cost[j,k];

        subj to OneJobPerMachine:
            alldiff {j in JOBS} MachineForJob[j];

        subj to CapacityOfMachine {k in MACHINES}:
            numberof k in ({j in JOBS} MachineForJob[j]) <= cap[k];

        subj to VisitOnce {j in BOATS}:
            isH[j] = 0 ==> alldiff {t in TIMES} H[j,t];

        subj to HostNever {j in BOATS}:
            isH[j] = 1 ==> forall {t in TIMES} H[j,t] = j;


- QP:

  .. code-block:: ampl

        -5 * (abs(x[1])-0.7)^2 + x[2]


- General non-linear:

  .. code-block:: ampl

        (x+1)*(x+2)*(x+3) + exp(y) == 8;


- Gurobi “general constraints” `and`, `or`, `max`, `min`, `abs` which
  have good MIP redefinitions, as well as `indicator (==>)`, are converted to
  MIP if finite bounds available, with tight or big-*M* constraints.
  This behavior can be changed by solver options **acc:abs** etc.
  to use Gurobi native constraints instead.

- Nonlinear “generals” are passed to Gurobi:

  `exp`, `log`, `sin`, `cos`, `tan`, `pow`, `pl`, `SOS1/2`.

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

and knowing that  `-15 <= x-y-z <= 30`, reformulate:

.. code-block:: ampl

        var t >=-15, <=30;
        subj to: t == x-y-z;
        subj to: log(x+2)<=y^2  or  t>=0;
