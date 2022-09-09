.. _getting-started:

Getting started
===============

.. _howto-create-new-driver:

Create a new driver
-------------------

For a driver setup of your choice (see :ref:`possible driver setups <driver-setups>`),
you can use existing drivers as templates. The process is detailed below.

Mock template driver 'visitor'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The easiest way to getting started developing a new solver driver using
`mp <https://github.com/ampl/mp>`_ is by
looking at the `visitor <https://github.com/ampl/mp/tree/develop/solvers/visitor>`_ mock
driver on branch ``develop``.

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



.. _driver-setups:

Solver driver setups
--------------------

A minimal setup
~~~~~~~~~~~~~~~

To write a minimal AMPL solver driver, it requires an
:ref:`NL file reader and .sol file writer <NL-SOL-files>`.

The recommended setup
~~~~~~~~~~~~~~~~~~~~~

The recommended driver structure is to use the
:ref:`Backend class hierarchy <backend-classes>`.
Creating such driver from a template is
:ref:`described here <howto-create-new-driver>`.

A legacy setup
~~~~~~~~~~~~~~

A legacy setup used :ref:`Solver classes <solver-classes>`.
Their support may be discontinued.
