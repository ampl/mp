How to...
=========



.. _howto-create-new-driver:

Create a new driver
-------------------


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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* First, clone the mp repository.
  Then either:

  #. Copy all the directory :file:`visitor` (or any other exiting driver files)
     into a new directory - and change its name.

  #. Rename all occurrences of the word "visitor".


  Or:

  #. Use the file createDriver.py, which does the two items above automatically


* Create and regularly use tests for essential functionality.
  See :ref:`howto-test`.

* Create a pull request.


.. _howto-test:

Test
-------

Ensure the major functionality of your driver is tested by providing end-to-end
tests, see below. You can add unit tests as well, see folder :file:`test`.


Run end-to-end tests
~~~~~~~~~~~~~~~~~~~~

* To **run all tests** for one or several ``solver``\ s, enter the command:

  .. code-block:: console

      python3 test/end2end/run.py solver [another_solver [...]]

  The ``solver``s and ``ampl`` executables are expected to be on the system path.

* Detailed results are saved into a CSV report file, see ``--reportstub``.

* To run a **subset of the test cases**, ``cd`` into the corresponding
  subfolder of ``test/end2end/cases``, or use the ``--dir`` or
  ``--nonrecursive`` options.


Add solver to the test library
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* To **add a new solver** to the test library, derive a corresponding class in
  ``Solver.py`` and list it in ``SolverCollection.py``.


Add new test cases
~~~~~~~~~~~~~~~~~~

The major functionality of a solver driver should be end-to-end tested.
Folder :file:`test/end2end/cases/categorized/fast` contains test cases
which can be run in a few seconds for this purpose.

* To **add new test cases**, add the model/data/AMPL script files and describe
  them in the local ``modellist.json`` having the following format. The top JSON
  object is an array of test cases. Each element is a dictionary with the
  following items, where non-compulsory items are italicized:

  * **"name": "<name>"**: <name> is the case name. Unless *"files"* is present
    (see below), the first word of <name> must coincide with the
    model / script stem. For example, test case using ``diet.mod`` could be
    called ``diet objno=5``.

  * **"tags": ["linear", "continuous"]**: tags specifying model type, the case
    is executed only if the tags are a subset of the solver's ones. Except the
    tag **"run"** which means the test case is an AMPL script.

  * *"files": ["diet.mod", "diet.dat"]*

  * *"objective": value*: expected objective value.

  * *"options": { "ANYSOLVER_options": "iisfind=1", "baron_options": "iisfind=12", "send_statuses": "0" }*.
    Option key ending with ``SOLVER_options`` is for any solver, except when
    a solver-specific key is present (like ``baron_options``.)

  * *"values": { "X[0].iis": "upp", ... }*. Expected values or expressions,
    in the form AMPL ``display`` command would accept.

    * For example, to check *logical expressions*, use if/then:

      .. code-block:: json

          "values": {
            "if color['Belgium'] != color['France'] then 1 else 0": 1,
            "solve_result_num": 0
          }


