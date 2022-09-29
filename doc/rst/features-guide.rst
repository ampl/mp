.. _features-guide:

Features Guide for MP-based AMPL Solvers
****************************************

.. highlight:: ampl

The mp framework defines standard solver *features* that solvers might support; 
these are usually characterized by a set of options used to control the feature,
sometimes suffixes to pass required data and results, and may change the behaviour
of the solution process.
Much of the biolerplate code is written already, so that the behaviour becomes
automatically standardized across all solvers.

This page presents the semantics of the most common solver features; for a development
reference refer to the `Visitor driver
<https://github.com/ampl/mp/tree/develop/solvers/visitor>`_ that shows how to
declare that the solver driver supports a feature and how to implement it.


General
=======


.. _warm-start:

Warm start
----------

Solution process can often benefit of a solution (a set of variable values) to start the algorithm. 
This is passed to supporting solver automatically if the option is activated and variables in AMPL
have a value assigned. Note that, for LP problems, also the dual values can be passed.

This option controls whether to use incoming primal (and dual, for LP) variable values in 
a warmstart.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``alg:start``
   * - **Applicability**
     - LP and MIP models
   * - **Input**
     - Variable values
   * - **Output**
     - None
   * - **Values**
     - Sum of:

       * **0** - No (default)
       * **1** - Yes (for LP: if there is no incoming alg:basis) (default)
       * **2** - Yes (for LP: ignoring the incoming alg:basis, if any)
   * - **Example**
     - Use :ref:`this model <nQueensModel>`

       Execute::

         option <solver>_options "alg:start=1";
         let n:=250; # increase the size of the model to have noticeable solution times

         solve;
         printf "Solution without warm start took %fs\n", _solve_time;

         # Now an optimal solution is already present, we pass it to the solver
         # which would use to start the solution process
         solve;
         printf "Solution with warm start took %fs\n", _solve_time;

       Output:

       .. code-block:: shell
 
            ...

            Solution without warm start took 2.89062s
            
            ...

            Solution with warm start took 0.671875s

.. _basisio:

Input and output basis
----------------------

A basis is a set of variable values representing a feasible and extreme solution.
Simplex solvers normally calculate this as part of the solution process, while
interior point methods must perform additional steps (crossover) to get it.
In a way similar to :ref:`warm start <warm-start>`, a basis can also be passed to the solver,
which will use it as starting point for searching for a solution.

This option controls whether to use or return a basis.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``alg:basis``
   * - **Applicability**
     - LP and MIP models
   * - **Input**
     - Suffix:

       * ``sstatus`` on variables 
   * - **Output**
     - Suffix:

       * ``sstatus`` on variables 
   * - **Values**
     - Sum of:

       * **0** - No
       * **1** - Use incoming basis (if provided)
       * **2** - Return final basis
       * **3** - Both (1 + 2, default)

   * - **Example**
     - Use :ref:`this model <multiObjectiveDiet>`

       Execute::

         option <solver>_options "alg:start=0 outlev=1"; # disable passing the solution

         solve;
         display Buy.sstatus; # display basis status

         solve; # second solve with take much less although a solution is not provided

       In the solver logs, we can see the expected behaviour:

       .. code-block:: shell

          x-Gurobi 9.5.2: optimal solution; objective 74.27382022
          3 simplex iterations
          Objective = total_cost['A&P']
          ampl: display Buy.sstatus;
          Buy.sstatus [*] :=
          BEEF  low
          CHK  upp
          FISH  low
          HAM  low
          MCH  low
          MTL  bas
          SPG  bas
          TUR  low;

          ... # second solve:          

          Solved in 0 iterations and 0.00 seconds (0.00 work units)
          

.. _feasibilityrelaxation:

Feasibility Relaxation
----------------------

The feasibility relaxation functionality enables the solver to find a feasible
solution even if the original model is unfeasible without explicitly adding
slack variables to the constraints.
In the feasibility relaxation problem, 

#. Each variable :math:`x` can violate its bounds (:math:`lb \leq x \leq ub`):
  
   * Violation of lower bound :math:`lbv = max(0, lb-x)`
   * Violation of upper bund :math:`ubv = max(0, x-ub)`

#. Each constraint body :math:`c` can violate its bounds also (:math:`c \leq rhs`)

   * Constraint violation :math:`rhsv = max(0, c-rhs)`

The objective then becomes to minimize some function of the
violations (e.g. the number of violations, or their sum - possibly weighted by some
penalty values).
The penalty values (used in some kinds of feasibility relaxation problmes) can be 
controlled with macro defaults (e.g. option ``alg:ubpen`` sets the penalty weight for 
all upper buonds violations, and its default values is 1) or, with more granularity,
on each entity via suffix values (e.g. variable suffix ``ubpen`` on variables, default
value 0). 
Penaly weights < 0 are treated as Infinity, allowing no violation.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``alg:feasrelax``
   * - **Applicability**
     - LP and MIP models
   * - **Input**
     -  * Options
  
          * ``alg:lbpen``: penalty for lower bound violations if suffix ``lbpen`` is not defined - default 1
          * ``alg:ubpen``: penalty for upper bound violations if suffix ``ubpen`` is not defined - default 1
          * ``alg:rhspen``: penalty for rhs violations if suffix ``rhspen`` is not defined - default 1

        * Suffixes

          * ``lbpen`` on variables - penalty for lower bound violations - default 0
          * ``ubpen`` on variables - penalty for upper bound violations - default 0
          * ``rhspen`` on constraints - penalty for rhs violations - default 0
   * - **Output**
     - None
   * - **Values**
     - Sum of:

       * **0** - No
       * **1** - Yes, minimizing the weighted sum of violations
       * **2** - Yes, minimizing the weighted sum of squared violations
       * **3** - Yes, minimizing the weighted count of violations
       * **3-6** - Same objective as 1-3, but also optimize the original objective, subject to the violation being minimized
   * - **Example**
     - Use :ref:`this model <infeasibleModel>`

       Solve the model changing the penalties to get different solutions::

          options <solver>_options "alg:feasrelax=1";
          option presolve 0;

          solve; display x,y;

       Gives:

        .. code-block:: bash

          x = 1
          y = 1
        
       Now we want to force all variable lower bounds to be respected::

          options <solver>_options "alg:feasrelax=1 alg:lbpen=-1";
          solve; display x,y;

       Gives, as expected x=5 (it had a lower bound of 5):

        .. code-block:: bash

            x = 5
            y = 1
        
       Single violations can be controlled via suffixes; in this case we
       want to control the constraints violations::

          options <solver>_options "alg:feasrelax=1 alg:lbpen=1"; # allow lower bounds to be violated again

          suffix rhspen IN;
          let C1.rhspen := 1; # normal weight
          let C2.rhspen := -1; # C2 can NOT be violated
          let C3.rhspen := 10; # We'd rather not violate C3
          solve;
          display C1.slack, C2.slack, C3.slack;

       Gives:

        .. code-block:: bash

          C1.slack = -3
          C2.slack = 4
          C3.slack = 0

       C1 is violated - which makes sense as we specified that C2 cannot be violated
       and we gave an higher avoidance weight to C2. If we want to violate C1 instead,
       we can::

          let C1.rhspen := 10; # We'd rather not violate C1
          let C2.rhspen := 1; # C2 can be violated
          let C3.rhspen := 1; # Normal weight
          solve;
          display C1.slack, C2.slack, C3.slack;

       Which gives, as expected:

        .. code-block:: bash

          C1.slack = 18
          C2.slack = 2
          C3.slack = -1


.. _multiplesolutions:

Multiple solutions
------------------

More often than not, optimization problems have more than one optimal solution; moreover, during the 
solution process, MIP solvers usually find sub-optimal solutions, which are normally discarded.
They can be however be kept, and in most cases there are solver-specific options to control how
the search for additional solutions is performed.

The main (and generic) options that controls the search are ``sol:stub`` amd ``sol:count``, which
control respecitvely the base-name for the files where additional solution will be stored and
if to count additional solutions and return them in the ``nsol`` problem suffix.
Specifying a stub name automatically enables the solutions count; found solutions are written to 
files [``solutionstub1.sol'``,  ... ``solutionstub<nsol>.sol``].


.. list-table::
   :header-rows: 0

   * - **Option**
     - ``sol:stub``
   * - **Applicability**
     - MIP models
   * - **Input**
     - None
   * - **Output**
     - Suffixes ``nsol`` and ``npool`` on problem
   * - **Values**
     - The name used as base file name for the alternative solutions

   * - **Example**
     - Use :ref:`this model <nQueensModel>`

       Execute::

         option <solver>_options "sol:stub=queentake";
         solve;
         printf "I have found %d solutions\n", Initial.nsol;

         printf "Displaying solution 1";
         solution queentake1.sol;
         display X;

         printf "Displaying solution 2";
         solution queentake2.sol;
         display X;

       Output:

       .. code-block:: shell
 
            x-Gurobi 9.5.2: sol:stub=queentake
            x-Gurobi 9.5.2: optimal solution; objective 10
            355 simplex iterations
            23 branching nodes

            suffix nsol OUT;
            suffix npool OUT;

            I have found 10 solutions
            Displaying solution 1
            x-Gurobi 9.5.2: Alternative solution; objective 10
            X [*,*]
            :    1   2   3   4   5   6   7   8   9  10    :=
            1    0   0   0   0   0   0   1   0   0   0
            2    1   0   0   0   0   0   0   0   0   0
            3    0   0   0   1   0   0   0   0   0   0
            4    0   0   0   0   0   1   0   0   0   0
            5    0   0   0   0   0   0   0   0   1   0
            6    0   0   1   0   0   0   0   0   0   0
            7    0   0   0   0   0   0   0   0   0   1
            8    0   0   0   0   0   0   0   1   0   0
            9    0   1   0   0   0   0   0   0   0   0
            10   0   0   0   0   1   0   0   0   0   0;

            Displaying solution 2
            x-Gurobi 9.5.2: Alternative solution; objective 10
            X [*,*]
            :    1   2   3   4   5   6   7   8   9  10    :=
            1    0   0   1   0   0   0   0   0   0   0
            2    0   0   0   0   0   1   0   0   0   0
            3    0   0   0   0   0   0   0   0   1   0
            4    1   0   0   0   0   0   0   0   0   0
            5    0   0   0   1   0   0   0   0   0   0
            6    0   0   0   0   0   0   1   0   0   0
            7    0   0   0   0   0   0   0   0   0   1
            8    0   0   0   0   0   0   0   1   0   0
            9    0   1   0   0   0   0   0   0   0   0
            10   0   0   0   0   1   0   0   0   0   0


.. _sensitivityAnalysis:

Sensitivity analysis
--------------------

It is often useful to know the ranges of variables and constraint bodies for which the optimal basis
remains optimal. Solvers supporting this feature return such ranges in suffixes after solving to optimum.
This option controls whether to calculate these values and return them in the suffixes listed below.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``alg:sens``
   * - **Applicability**
     - LP and MIP models
   * - **Input**
     - None
   * - **Output**
     - Suffix:
       * ``sensobjlo`` variables, smallest objective coefficient
       * ``sensobjhi`` variables, greatest objective coefficient
       * ``senslblo`` variables, smallest variable lower bound
       * ``senslbhi`` variables, greatest variable lower bound
       * ``sensublo`` variables, smallest variable upper bound
       * ``sensubhi`` variables, greatest variable upper bound
   * - **Values**
     - Sum of:

       * **0** - No (default)
       * **1** - Yes
       
   * - **Example**
     - Use :ref:`multiObjectiveDiet`

       Execute::

          options <solver>_options "alg:sens=1"; 
          solve;

       Then the ranges for variables and constraints can be examined: :

          display Buy, Buy.sstatus, Buy.sensublo, Buy.sensubhi;

       Which gives:

       .. code-block:: bash

          display Buy.sensublo, Buy.sensubhi, Buy.sstatus, Buy;
          :    Buy.sensublo  Buy.sensubhi Buy.sstatus     Buy       :=
          BEEF    2           1e+100        low          2
          CHK     9.12987         10.9792   upp         10
          FISH    2           1e+100        low          2
          HAM     2           1e+100        low          2
          MCH     2           1e+100        low          2
          MTL     6.23596     1e+100        bas          6.23596
          SPG     5.25843     1e+100        bas          5.25843
          TUR     2           1e+100        low          2;          

.. _kappa:

Kappa
-----

Kappa is the condition number for the current LP basis matrix. 
It is a measure of the stability of the current solution :math:`Ax=b`
measuring the rate at which the solution :math:`x` will change with respect to a 
change in :math:`b`. 
It is only available for basic solutions, therefore it is not available for barrier method
if crossover is not applied.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``alg:kappa``
   * - **Applicability**
     - LP and MIP models with optimal basis
   * - **Input**
     - None
   * - **Output**
     - Additional text in ``solve_message`` and suffix:

       * ``kappa`` on objective and problem
   * - **Values**
     - Sum of:

       * **0** - No
       * **1** - Report kappa in solve_message
       * **2** - Return kappa in the solver-defined suffix ``kappa``
   * - **Example**
     - Use :ref:`multiObjectiveDiet`

       Solve the model and report kappa:

          options <solver>_options "alg:kappa=3"; 
          solve;

          display Initial.kappa, total_number.kappa;
 
       Gives:

       .. code-block:: bash

        x-Gurobi 9.5.2: optimal solution; objective 30.92537313
        kappa value: 53.5399
        5 simplex iterations

        suffix kappa OUT;

        Initial.kappa = 53.5399
        total_number.kappa = 53.5399


.. _unboundedRays:

Unbounded rays
--------------

When a model is unbounded, a vector :math:`r` (unbounded ray)  can be found such that
when added to any feasible solution :math:`x`, the resulting vector is a feasible solution
with an improved objective value.
When a model is infeasible, the dual solution is unbounded and the same as above can be applied
to constraints.
This option controls whether to return suffix ``unbdd`` if the objective is unbounded
or suffix ``dunbdd`` if the constraints are infeasible.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``alg:rays``
   * - **Applicability**
     - Unbounded/unfeasible LP and MIP linear models
   * - **Input**
     - None
   * - **Output**
     - Suffixes:

       * ``unbdd`` on variables if the problem is unbounded
       * ``dunbdd`` on constraints if the problem is infeasible
   * - **Values**
     - Sum of:

       * **0** - Do not calculate or return unbounded rays
       * **1** - Return only ``unbdd``
       * **2** - Return only ``dunbdd``
       * **3** - Return both (default)
   * - **Example**
     - Use :ref:`this model <iisModel>`

       Solve the (infeasible) model and report the ``dunbdd`` rays::

          options <solver>_options "alg:rays=3"; # it is default already
          solve;
          display c1.dunbdd, c2.dunbdd, c3.dunbdd;

       Output:

       .. code-block:: shell

          x-Gurobi 9.5.2: alg:rays=3
          x-Gurobi 9.5.2: infeasible problem

          suffix dunbdd OUT;

          c1.dunbdd = 1
          c2.dunbdd = 0
          c3.dunbdd = 0


.. _multipleObjectives:

Multiple objectives
-------------------

Many real world problems have multiple objectives; often this scenario is tackled by blending all the objectives
by linear combination when formulating the model, or by minimizing each unwanted objective deviations from a pre-specified
goal.
Many solvers can facilitate the formulation; the available functionalities are solver-specific; at MP level
they accessible via the main option ``obj:multi``. Consult the solver documentation for the functionalities available
on your solver.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``obj:multi``
   * - **Applicability**
     - LP and MIP models
   * - **Input**
     - None
   * - **Output**
     - All objectives are reported in ``solve_message``

   * - **Values**
     - Values:

       * **0** - No (default)
       * **1** - Yes
   * - **Example**
     - Use :ref:`multiObjectiveDiet`

       Execute::

          options <solver>_options "obj:multi=1"; 
          solve;

       Output:

       .. code-block:: shell

          x-Gurobi 9.5.1: obj:multi=1
          x-Gurobi 9.5.1: optimal solution; objective 74.27382022
          Individual objective values:
            _sobj[1] = 74.27382022
            _sobj[2] = 75.01966292
            _sobj[3] = 79.59719101
            _sobj[4] = 31.49438202


   
Irreducible Inconsistent Set (IIS)
----------------------------------

Given an infeasible model, it is useful to know where the infeasibility comes from, that is, which
bounds and/or constraints are incompatible.
An IIS is a subset of contraint and variables that is still infeasible and where if a member is removed
the subsystem becaomes feasible. Note that an infeasible model can have more than one IIS.

This options controls whether to perform the additional computational steps required to find an IIS.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``alg:iisfind``, ``iisfind``,  ``iis``
   * - **Applicability**
     - LP and MIP models
   * - **Input**
     - None
   * - **Output**
     - Suffix ``iis`` on variables and constraints

   * - **Values**
     - Values:

       * **0** - No (default)
       * **1** - Yes
   * - **Example**
     - Use :ref:`this model <iisModel>`

       Execute::

        option <solver>_options "iisfind=1";
        option presolve 0; # else the model would be oversimplified
        solve;
        display c1.iis, c2.iis, c3.iis;

       Output:

       .. code-block:: shell

          x-Gurobi 9.5.2: alg:iisfind=1
          x-Gurobi 9.5.2: infeasible problem
          2 simplex iterations

          suffix iis symbolic OUT;
          suffix dunbdd OUT;
  
          c1.iis = mem
          c2.iis = mem
          c3.iis = non


MIP only
========


.. _returnMIPgap:

Return MIP gap
--------------

The MIP gap describes how far the reported solution for the MIP model is from the
best bound found. It is defined as:

:math:`absgap = | objective - bestbound |`

:math:`relmipgap = absgap / | objective |``

It gives a measure of how far the current solution is from the
theoretical optimum.
The AMPL option controls whether to return mipgap suffixes or include mipgap values 
in the solve_message. Returned suffix values are ``+Infinity`` if no integer-feasible 
solution has been found, in which case no mipgap values are reported in the solve_message.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``mip:return_gap``
   * - **Applicability**
     - MIP models
   * - **Input**
     - None
   * - **Output**
     - Additional text in ``solve_message`` and suffixes:

       * ``relmipgap`` on objective and problem
       * ``absmipgap`` on objective and problem
   * - **Values**
     - Sum of:

       * **0** - Do not report gaps
       * **1** - Return .relmipgap suffix (relative to ``|obj|``)
       * **2** - Return .absmipgap suffix (absolute mipgap)
       * **4** - Suppress mipgap values in solve_message
   * - **Example**
     - Use :ref:`this model <nQueensModel>`

       Execute::

         # 3 = return both absolute and relative MIP gap, and also report them
         # in the solve_message
         option <solver>_options "return_mipgap=3";
         solve;

         display max_queens.relmipgap, max_queens.absmipgap;
         display Initial.relmipgap, Initial.absmipgap;

       Output:

       .. code-block:: bash

          x-Gurobi 9.5.2: mip:return_gap=3
          x-Gurobi 9.5.2: optimal solution; objective 10

          suffix absmipgap OUT;
  
          max_queens.relmipgap = 0
          max_queens.absmipgap = 0
          Initial.relmipgap = 0
          Initial.absmipgap = 0



.. _returnBestBound:

Return best dual bound
----------------------

The best dual bound (on the objective) represents what is the currently 
proven best value that the objective value can assume. Usually solvers terminate
when the current solution is close enough to the best bound.

This option controls whether to return suffix .bestbound for the best known MIP dual
bound on the objective value. The returned value is -Infinity for minimization
problems and +Infinity for maximization problems if there are no integer 
variables or if a dual bound is not available.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``mip:bestbound``
   * - **Applicability**
     - MIP models
   * - **Input**
     - None
   * - **Output**
     - Suffix:

       * ``bestbound`` on objective
   * - **Values**
     - Sum of:

       * **0** - No (default)
       * **1** - Yes
   * - **Example**
     - Use :ref:`this model <nQueensModel>`

       Execute::

            option <solver>_options "mip:bestbound=3";
            solve;
            display max_queens.bestbound, Initial.bestbound;

       Output:

       .. code-block:: shell

            x-Gurobi 9.5.2: mip:bestbound=1
            x-Gurobi 9.5.2: optimal solution; objective 10
            
            max_queens.bestbound = 10
            Initial.bestbound = 10

.. _lazyConstraints:

Lazy constraints and user cuts
------------------------------

The solution process of a MIP model can be helped by further specifying its structure.
Specifically, constraints can be marked as ``lazy`` or as ``user cuts``.
Such constraints are not included initially, then:

``lazy cosntraints`` are pulled in when a feasible solution is found; if the solution violates
them, it is cut off. They are integral part of the model, as they can cut off integer-feasible 
solutions.

``user cuts`` are pulled in also, but they can only cut off relaxation solutions; they are 
consider redundant in terms of specifying integer feasibility.

This option controls whether to recognize the suffx ``lazy`` on constraints, which should then be 
positive to denote a lazy constraint and negative to mark a contraint as a user cut.


.. list-table::
   :header-rows: 0

   * - **Option**
     - ``mip:lazy``
   * - **Applicability**
     - MIP models
   * - **Input**
     - Suffix:

       * ``lazy`` on constraints (>0 for lazy constraint, <0 for user cuts) 
   * - **Output**
     - None
  
   * - **Values**
     - Sum of:

       * **0** - No
       * **1** - Accept >0 values to denote lazy constraints
       * **2** - Accept <0 values to denote user cuts
       * **3** - Accept both (default)
   * - **Example**
     - Following ampl model

       .. code-block:: ampl
 
            TODO Model

       end of code block

       .. code-block:: ampl

            suffix priority IN;
            let c1.priority := 1;  # lazy constraint
            let c2.priority := -1; # user cut, must be redundant as it might never be pulled in
            solve;

            # TODO SHOW OUTPUT



.. _varPriorities:

Variable priorities
-------------------

Solution of MIP models via branch and bound can often be helped by providing
preferences on which variables to branch on. Those can be specified in AMPL via the suffix 
``priority``.

This option controls whether to read those values and use them in the solution process.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``mip:priorities``
   * - **Applicability**
     - MIP models
   * - **Input**
     - Suffix:

       * ``priority`` on variables 
   * - **Output**
     - None
   * - **Values**
     - Values:

       * **0** - Ignore priorities
       * **1** - Read priorities (default)
   * - **Example**
     - Following AMPL model

       .. code-block:: ampl
 
            TODO Model

       end of code block

       .. code-block:: ampl

            let x.priority := 1;
            let y.priority := 5;
            solve;

            # TODO SHOW OUTPUT




.. _fixedModel:

Fixed model (return basis for MIP)
----------------------------------

At the end of the solution process for a MIP model, a continuous relaxation of the model
with all the integer variables fixed at their integer-optimum value. Some continuous variables
can also be fixed to satisfy SOS or general constraints.
The model can therefore be solved without these types of restrictions to calculate a basis,
dual values or sensitivity information that wouldn't normally be available for MIP problems.

This option controls if to generate and solve the fixed model after solving the integer problem.

.. list-table::
   :header-rows: 0

   * - **Option**
     - ``mip:basis``
   * - **Applicability**
     - MIP models
   * - **Input**
     - None

   * - **Output**
     - Suffixes:
  
       * ``dual`` on variables
       * ``ssttus`` on variables
       * See :ref:`sensitivityanalysis` if requested
  
   * - **Values**
     - Values

       * **0** - No (default)
       * **1** - Yes
   * - **Example**
     - Following ampl model

       .. code-block:: ampl
 
            TODO Model

       end of code block

       .. code-block:: ampl

            option <solver>_options "mip:basis=1";

            # TODO SHOW OUTPUT


* Round


Reference models
================
.. toctree::
   :maxdepth: 2

    Reference models <models>

