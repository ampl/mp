.. _features-guide:

Features Guide for MP-based AMPL Solvers
****************************************

The mp framework defines standard solver *features* that solvers might support; 
these are usually characterized by a set of options used to control the feature,
sometimes suffixes to pass required data and results, and may change the behaviour
of the solution process.
Much of the biolerplate code is written alrady, so that the behaviour becomes 
automatically standardized across all solvers.

This page presents the semantics of the most common solver features; for a development
reference refer to the `Visitor driver
<https://github.com/ampl/mp/tree/develop/solvers/visitor>`_. that shows how to
declare that the solver driver supports a feature and how to implement it.


General
=======

* Kappa
* :ref:`feasibiliyRelaxation`
* Multiple objectives
* :ref:`multipleSolutions`
* :ref:`returnBestBound`
* :ref:`basisio`
* Return rays
* IIS (return and/or force)
* Sensitivity analysis

MIP only
========

* :ref:`returnMIPgap`
* :ref:`returnBestBound`
* Lazy constraints / user cuts
* Variable priorities // TODO Plural??
* Fixed model (return basis for MIP)
* Round

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
     - Following ampl model

       .. code-block:: ampl
 
            indented because i amp 
            a code block

       end of code block

       .. code-block:: ampl

            # 3 = return both absolute and relative MIP gap, and also report them
            # in the solve_message
            option <solver>_options "return_mipgap=3";
            solve;

            display obj.relmipgap, obj.absmipgap;
            display Initial.relmipgap, Initial.absmipgap;
            display solve_message;

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
     - Following ampl model

       .. code-block:: ampl
 
            TODO model

       TODO then

       .. code-block:: ampl

            option <solver>_options "mip:bestbound=3";
            solve;

            display obj.bestbound;


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
     - Following ampl model

       .. code-block:: ampl
 
            indented because i amp 
            a code block

       end of code block

       .. code-block:: ampl

            option <solver>_options "alg:start=1";
            let var x := // TODO set to actual solution
            solve;




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
     - Following ampl model

       .. code-block:: ampl
 
            TODO Model

       end of code block

       .. code-block:: ampl

            option <solver>_options "alg:basis=3"; # already set by default
            display x.sstatus // TODO
            solve;

            # TODO SHOW OUTPUT


.. _feasibiliyrelaxation:

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
     - Following ampl model

       .. code-block:: ampl
 
          # TODO Infeasible model

       Solve the model changing the penalties to get different solutions: 

       .. code-block:: ampl

          # No lower bound can be violated
          options <solver>_options "lbpen=-1";
          suffix rhspen IN;

          let C1.rhspen := 1; # normal weight
          let C2.rhspen := -1; # C2 can NOT be violated
          let C3.rhspen := 10; # We'd rather not violate C3
          let C4.rhspen := 0; # We don't care if we violate C4

          solve;

          display C1.slack, C2.slack, C3.slack, C4.slack;

          let C2.rhspen := 1; # C2 can be violated
          let C3.rhspen := 10; # We'd rather not violate C3
          let C4.rhspen := 0; # We don't care if we violate C4

          display C1.slack, C2.slack, C3.slack, C4.slack;


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