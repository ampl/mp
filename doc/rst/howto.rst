.. _howto:

HOWTO hook your solver to AMPL MP
=================================

This page describes the following steps to create a new driver:

.. contents::
   :depth: 1
   :local:
   :backlinks: none


A reference of the API is :ref:`here <components>`.


.. _driver-setups:

Solver driver setups
--------------------

.. _driver-minimal-setup:

A minimal setup
~~~~~~~~~~~~~~~

To write a minimal AMPL solver driver, it requires an
:ref:`NL file reader and .sol file writer <NL-SOL-files>`.
Examples:
`LocalSolver <https://github.com/ampl/mp/tree/develop/solvers/localsolver>`_,
`SCIP 8.0 <https://scipopt.org/>`_.

.. _driver-recommended-setup:

Recommended setup
~~~~~~~~~~~~~~~~~~~~~

The recommended driver structure is to use the
:ref:`Backend class hierarchy with FlatModelAPI <recommended-driver-logic>`.
Examples: :ref:`flat-solvers`.


A legacy setup
~~~~~~~~~~~~~~

A legacy setup used :ref:`Solver classes <solver-classes>`.
Their support may be discontinued. Examples include
`Ilogcp <https://github.com/ampl/mp/tree/develop/solvers/ilogcp>`_
and
`Gecode <https://github.com/ampl/mp/tree/develop/solvers/gecode>`_.


.. _howto-create-new-driver-from-template:

Clone a new driver from a template
-----------------------------------

For a driver setup of your choice (see :ref:`possible driver setups <driver-setups>`),
you can use existing drivers as templates. The process is detailed below.

Mock template driver 'visitor'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The easiest way to getting started developing a new solver driver using
`mp <https://github.com/ampl/mp>`_ is by
looking at the `visitor <https://github.com/ampl/mp/tree/develop/solvers/visitor>`_ mock
driver on branch ``develop``. This template uses
:ref:`Backend class hierarchy <backend-classes>` with :ref:`flat-model-api`.

To build it, you can configure the build system of your choice and specify
the cmake variable `BUILD` appropriately::

  mkdir build
  cd build
  cmake .. -DBUILD=visitor
  make

For faster recompilation, install ``ccache`` and
add the following CMake flags::

  -DBUILD_TESTS=off -DBUILD_EXAMPLES=off -DBUILD_DOC=off
  -DCMAKE_BUILD_TYPE=Debug                     ## Linux/Unix way to set debug mode
  -DUSE_SANITIZERS=on                          ## Linux/Unix way to use code sanitizers
                                               ## (slow, for checking only)

Once built, executing::

  ./visitor modelfilename.nl

will execute the mock driver, which will simply visit the model represented
in the nl file.
The visitor source code can be used as a template to create a new driver,
as described in the section below.


Copying a driver template
~~~~~~~~~~~~~~~~~~~~~~~~~

* First, clone the mp repository.
  Then either:

  #. Copy all the directory :file:`solvers/visitor` (or any other exiting driver files)
     into a new directory - and change its name.

  #. Rename all occurrences of the word "visitor".


  or:

  #. Use the script :file:`solvers/createDriver.py`, which does the two items above
     automatically. The script expects a source driver and a new driver name. So,
     to create a new driver named ``brandNewAMPLSolver`` based on ``visitor``, execute::

        python3 createDriver.py visitor brandNewAMPLSolver


* Add the new target in :file:`solvers/CMakeLists.txt`.

* Implement :ref:`driver features <implement-basic-driver-features>`.

* Create a pull request.


.. _implement-basic-driver-features:

Get a basic version up & running
--------------------------------

To implement a barebone driver, several methods need to be specialized.
The description assumes a
:ref:`cloned visitor driver template <howto-create-new-driver-from-template>`,
which uses the :ref:`recommended driver setup <driver-recommended-setup>`.


.. _backend-vs-modelapi:

...Backend vs ...ModelAPI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As detailed in :ref:`recommended-driver-logic`, solver API is addressed by two wrapper objects,
:ref:`Custom Backend <backend-classes>` and :ref:`flat-model-api`.
Their names are normally <yourSolver>Backend and <yourSolver>ModelAPI, respectively.
While they perform different tasks, they have some common core, for example to keep the
underlying solver API pointer. This is managed by the common ancestor <yourSolver>Common:

- SolverCommon --> SolverBackend
- SolverCommon --> SolverModelAPI

(--> means inheritance). Thus, there are two objects of SolverCommon keeping the same
underlying solver API pointers (typically the environment and model pointers).
In fact, the information duplicated  between the two objects is stored in the extra class
SolverCommonInfo. Example:

.. code-block:: c++

   /// Information shared by both
   /// `HighsBackend` and `HighsModelAPI`
   struct HighsCommonInfo {
     void* lp() const { return lp_; }
     void set_lp(void* lp) { lp_ = lp; }
   private:
     void*      lp_ = nullptr;
   };


.. _basic-spec-model-api:

Basic specialization of ModelAPI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :ref:`ModelAPI wrapper <backend-vs-modelapi>` requires a minimal specialization
to be able to accept MILP models. Note that linear models still cover much of the
general modeling
capabilities highlighted in the :ref:`Modeling Guide <modeling-guide>`
(via automatic transformations). <yourSolver>ModelAPI should declare and implement
the following:

.. code-block:: c++

   /// Called before problem input.
   /// Model info can be used to preallocate memory.
   void InitProblemModificationPhase(const FlatModelInfo*);
   /// After
   void FinishProblemModificationPhase();

   void AddVariables(const VarArrayDef& );
   void SetLinearObjective( int iobj, const LinearObjective& lo );
   /// Whether accepting quadratic objectives:
   /// 0 - no, 1 - convex, 2 - nonconvex
   static int AcceptsQuadObj() { return 1; }
   void SetQuadraticObjective(int iobj, const QuadraticObjective& qo);

   //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
   USE_BASE_CONSTRAINT_HANDLERS(BaseModelAPI)

   ACCEPT_CONSTRAINT(LinConRange, Recommended, CG_Linear)
   void AddConstraint(const LinConRange& lc);
   ACCEPT_CONSTRAINT(LinConLE, Recommended, CG_Linear)
   void AddConstraint(const LinConLE& lc);
   ACCEPT_CONSTRAINT(LinConEQ, Recommended, CG_Linear)
   void AddConstraint(const LinConEQ& lc);
   ACCEPT_CONSTRAINT(LinConGE, Recommended, CG_Linear)
   void AddConstraint(const LinConGE& lc);


Note that after `InitProblemModificationPhase()` has been called,
you can access the `mp::FlatModelInfo` object with the inherited method
`GetFlatModelInfo()`.
For more advanced modeling, see :ref:`configure-automatic-model-conversions`.


Basic specialization of the Backend
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :ref:`Backend <backend-vs-modelapi>` requires a minimal specialization to
enable some very common :ref:`AMPL driver logic <features-guide>`.
<yourSolver>Backend should declare and
implement the following:

.. code-block:: c++

   /// Prefix used for the <prefix>_options environment variable
   static const char* GetAMPLSolverName() { return "highs"; }

   /// Solver name displayed in messages
   static const char* GetSolverName() { return "HiGHS"; }
   /// AMPL driver name displayed in messages
   static const char* GetAMPLSolverLongName() { return "AMPL-HiGHS"; }
   /// Version displayed with -v
   std::string GetSolverVersion();

   /// Init custom driver options, such as outlev, writeprob
   void InitCustomOptions() override;
   /// Chance for the Backend to init solver environment, etc
   void InitOptionParsing() override;
   /// Chance to consider options immediately (open cloud, etc)
   void FinishOptionParsing() override;

   /// Note the interrupt notifier
   void SetInterrupter(mp::Interrupter* inter) override;

   /// Solve, no model modification any more (such as feasrelax).
   /// Can report intermediate results via HandleFeasibleSolution() during this,
   /// otherwise/finally via ReportResults()
   void Solve() override;

   /// Report final solution
   void ReportHIGHSResults() override;

   /// Values of all objectives
   ArrayRef<double> GetObjectiveValues() override;
   /// Primal solution values. Empty if not available
   ArrayRef<double> PrimalSolution() override;
   /// Dual solution. Empty if not available
   pre::ValueMapDbl DualSolution() override;

   /// Solution attributes
   double NodeCount() const;
   double SimplexIterations() const;
   int BarrierIterations() const;

   /// Convert solution/solver status to code+string
   std::pair<int, std::string> ConvertHIGHSStatus();
   /// Add custom messages
   void AddHIGHSMessages();

For other common and custom features, see :ref:`implement-standard-features`.


.. _configure-automatic-model-conversions:

Configure automatic model reformulations
------------------------------------------

This section describes configuration of the
:ref:`automatic model reformulations <modeling-guide>`
provided by the AMPL MP library, as well as adding new reformulations.

Configure automatic reformulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The MP library can reformulate most of AMPL's
:ref:`nonlinear and logical expressions <modeling-guide>`
to basic MILP or MIQP constructs. If your solver does not natively handle
an expression, you don't have to code anything. But assuming the indicator
constraints are supported, the following code needs to be added:

.. code-block:: c++

  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinLE& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinEQ, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinEQ& mc);
  ACCEPT_CONSTRAINT(IndicatorConstraintLinGE, Recommended, CG_General)
  void AddConstraint(const IndicatorConstraintLinGE& mc);

If you want a big-M linearization to be attempted first, replace `Recommended` by
`AcceptedButNotRecommended`.
To see the list of supported constraints, which largely correspond to
:ref:`AMPL modeling expressions <modeling-guide>`, see `constr_std.h`,
or run an existing driver with `-c`.
For explanation of constraint groups, see :ref:`value-presolver`.

Specifically for quadratic constraints
(quadratic objectives were discussed in :ref:`basic-spec-model-api`),
implement

.. code-block:: c++

   /// Ask if the solver accepts non-convex quadratic constraints
   static constexpr bool AcceptsNonconvexQC() { return false; }

   /// QuadConRange is optional.
   ACCEPT_CONSTRAINT(QuadConRange, Recommended, CG_Quadratic)
   void AddConstraint(const QuadConRange& qc);

   /// If using quadratics,
   /// QuadCon(LE/EQ/GE) should have 'Recommended'
   /// and have an implementation.
   ACCEPT_CONSTRAINT(QuadConLE, Recommended, CG_Quadratic)
   void AddConstraint(const QuadConLE& qc);
   ACCEPT_CONSTRAINT(QuadConEQ, Recommended, CG_Quadratic)
   void AddConstraint(const QuadConEQ& qc);
   ACCEPT_CONSTRAINT(QuadConGE, Recommended, CG_Quadratic)
   void AddConstraint(const QuadConGE& qc);


Convex quadratic solvers can be used to solve nonconvex problems
via piecewise-linear approximation of quadratics. To force the approximation,
set options *cvt:quadobj=0 cvt:quadcon=0*.


.. _implement-new-model-conversions:

Add new model reformulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section describes how to add new model reformulations
in the :ref:`recommended driver setup <driver-recommended-setup>`.

An overview of the reformulation process is provided in
:ref:`mm-and-reformulations`.


Derive a custom FlatConverter (NOT recommeded)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is recommended not to derive a solver-specific FlatConverter class.
Instead, add your conversions to a standard class and make it
optional (see existing options, such as *cvt:quadcon*), making it available
to other solvers.

Only if it's very specific,
derive a custom class from `mp::FlatConverter` or `mp::MIPFlatConverter`
and use it in ``Create<YourSolver>ModelMgr`` (its default implementation
is in `<yourSolver>modelapi.cc`).


Add a converter for a constraint
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To add a new converter, derive a new class from
`mp::BasicItemConverter` or `mp::BasicFuncConstrCvt` and follow the pattern
of existing conversions, in particular the
``INSTALL_ITEM_CONVERTER`` macro.


Add a new constraint type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To add a new constraint type, follow the definitions and installation pattern
from the ``STORE_CONSTRAINT_TYPE__WITH_MAP`` and
``STORE_CONSTRAINT_TYPE__NO_MAP`` in `mp::FlatConverter`.
In particular, overload ``PreprocessConstraint`` and
``PropagateResult`` for the new type in your custom class.


.. _implement-standard-features:

Implement standard & custom features
----------------------------------------

This section describes implementation of the
:ref:`optional standard driver features <features-guide>`,
as well as solver-specific features.
Some standard features are very common, such as BASIS,
others not, such as FIX_MODEL,
and don't have to be implemented unless the solver directly supports them.
The workflow relies on the
:ref:`Backend class hierarchy <backend-classes>`.


General standard features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Output level
^^^^^^^^^^^^

To implement the :ref:`standard behaviour of option outlev <outlev>`,
do the following:

1. Add solver option *outlev*. Its values can be solver-specific but ideally
   0 means silent and values above 0 mean some verbosity. Example code:

   .. code-block:: c++

      AddSolverOption("tech:outlev outlev",
        "0*/1: Whether to write mosek log lines to stdout.",
        MSK_IPAR_LOG, 0, 1);

2. In method `OpenSolver()` set verbosity level to silent, before the options
   are processed.

3. In `FinishOptionParsing()` call the inherited method `set_verbose_mode(v)`
   with `v==true` iff *outlev>0*.


Sensitivity analysis
^^^^^^^^^^^^^^^^^^^^

To implement the :ref:`standard behavior of option sens <sensitivityAnalysis>`,
do the following:

1. In your `Backend` class, declare:

   .. code-block:: c++

      ALLOW_STD_FEATURE(SENSITIVITY_ANALYSIS, true)

2. For derivatives of `mp::FlatBackend` you can override `GetSensRangesPresolved()`
   which automatically :ref:`postsolves <>` the sensitivity information:

   .. code-block:: c++

      SensRangesPresolved GetSensRangesPresolved() override;

   Currently this requires the vectors *con(lb/ub)(lo/hi)* to be populated for all
   linear constraints, including *LinCon(LE/EQ/GE)*. See the MOSEK driver for
   an example.

3. Alternatively, override `GetSensRanges()`:

   .. code-block:: c++

      SensRanges GetSensRanges() override;

   and implement it so that it returns postsolved information. See the Gurobi driver
   for an example.


MIP-only standard features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Fixed model (return basis for MIP)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To implement the
:ref:`standard behavior of the (probably exotic) option mip:basis / fixmodel <fixedModel>`,
do the following:

1.  In your `Backend` class, declare:

   .. code-block:: c++

      ALLOW_STD_FEATURE( FIX_MODEL, true )

2. Check method `need_fixed_MIP()` which returns true of user wants the fixed MIP
   information. In this case, your implementation should fix all non-continuous
   variables and variables from SOS / piecewise-linear constraints
   to their optimal values and solve the resulting LP; subsequent calls
   to `GetBasis()`, as well as dual solution and sensitivity information should
   correspond to that LP solution.


.. _implement-custom-features:

Custom features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _implement-pre-postsolving:

Pre- and postsolving of solutions and suffixes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For API details, see :ref:`value-presolver`.


