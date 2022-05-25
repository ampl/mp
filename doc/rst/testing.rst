
.. _howto-test:

Testing
=======

Ensure the major functionality of your driver is tested by providing and using
end-to-end tests. Unit tests can be used as well.


End-to-end tests
----------------

End-to-end tests have proven most important.

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


Unit tests
----------

You can employ unit tests as well, see folder :file:`test`.
