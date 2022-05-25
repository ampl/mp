How to...
=========

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

  -DBUILD_TESTING=off -DBUILD_DOC=off
  -DCMAKE_BUILD_TYPE=Debug                     ## Linux/Unix way to set debug mode

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

  #. Copy all the directory :file:`visitor` (or any other exiting driver files)
     into a new directory - and change its name.

  #. Rename all occurrences of the word "visitor".


  Or:

  #. Use the file createDriver.py, which does the two items above automatically

* Add the new target in :file:`solvers/CMakeLists.txt`.

* Create and regularly use tests for essential functionality.
  See :ref:`howto-test`.

* Create a pull request.



.. _driver-setups:

Solver driver setups
--------------------

A minimal setup
~~~~~~~~~~~~~~~

To write a minimal AMPL solver driver, it requires an NL file reader
and .sol file writer.

The recommended setup
~~~~~~~~~~~~~~~~~~~~~

A legacy setup
~~~~~~~~~~~~~~

