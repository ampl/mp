.. _howto:

HOWTO develop a new driver
==========================

This page describes possible driver setups, adding standard features
and configuring model conversions, and using existing driver code as
a template.


.. _driver-setups:

Solver driver setups
--------------------

.. _driver-minimal-setup:

A minimal setup
~~~~~~~~~~~~~~~

To write a minimal AMPL solver driver, it requires an
:ref:`NL file reader and .sol file writer <NL-SOL-files>`.

.. _driver-recommended-setup:

The recommended setup
~~~~~~~~~~~~~~~~~~~~~

The recommended driver structure is to use the
:ref:`Backend class hierarchy <backend-classes>` with :ref:`flat-model-api`.

A legacy setup
~~~~~~~~~~~~~~

A legacy setup used :ref:`Solver classes <solver-classes>`.
Their support may be discontinued.


.. _implement-standard-features:

Implementing standard features
----------------------------------------

This section describes implementation of the
:ref:`standard driver features <features-guide>`.
This relies on the
:ref:`Backend class hierarchy <backend-classes>`.


General features
~~~~~~~~~~~~~~~~

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


MIP-only features
~~~~~~~~~~~~~~~~~


Fixed model (return basis for MIP)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To implement the :ref:`standard behavior of option mip:basis / fixmodel <fixedModel>`,
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


.. _implement-pre-postsolving::

Pre- and postsolving of solutions and suffixes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For API details, see :ref:`value-presolver`.


.. _configure-automatic-model-conversions:

Configuring automatic model conversions
---------------------------------------

This section describes configuration of the
:ref:`automatic model conversions <modeling-guide>`
provided by the AMPL MP library.


.. _implement-new-model-conversions:

Implementing new model conversions
----------------------------------

This section describes how to add new model conversions
to the `mp::FlatConverter` and related / derived classes.


.. _howto-create-new-driver-from-template:

Create a new driver from a template
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

* Create and regularly use tests for essential functionality.
  See :ref:`howto-test`.

* Create a pull request.


