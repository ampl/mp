
.. _howto-test:

Testing
=======

Ensure the major functionality of your driver is tested by :ref:`running
end-to-end tests <run_end2end_tests>`. :ref:`Add your own test cases
<add_new_test_cases>` for essential functions.
:ref:`Unit tests <unit_tests>` can be used as well.


End-to-end tests
----------------

End-to-end tests have proven most important.


.. _run_end2end_tests:

Run end-to-end tests
~~~~~~~~~~~~~~~~~~~~

* To **run all tests** for one or several ``solver``\ s, enter the command:

  .. code-block:: console

      python3 test/end2end/run.py solver [another_solver [...]]

  The ``solver`` and ``ampl`` executables are expected to be on the system path.

* Detailed results are saved into a CSV report file, see ``--reportstub``.

* To run a **subset of the test cases**, ``cd`` into the corresponding
  subfolder of ``test/end2end/cases``, or use the ``--dir`` or
  ``--nonrecursive`` options.


Add solver to the test library
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* To **add a new solver** to the test library, derive a corresponding class in
  ``Solver.py`` and list it in ``SolverCollection.py``.

  * The *stags* parameter of the solver definition base class constructor
    in ``Solver.py``
    should contain the tags describing available driver features.


.. _add_new_test_cases:

Add new test cases
~~~~~~~~~~~~~~~~~~

The major functionality of a solver driver should be end-to-end tested.
Folder :file:`test/end2end/cases/categorized/fast` contains test cases
which can be run in a few seconds for this purpose, which should be done
frequently.

To **add new test cases**, add the model/data/AMPL script files in
a subfolder of :file:`test/end2end/cases/categorized/` and describe
them in the local ``modellist.json`` having the following format.
The top JSON
object is an array of test cases. Each element is a dictionary with the
following items, where non-compulsory items are italicized:

* **"name": "<name>"**: case name. Unless *"files"* is present
  (see below), the first word of <name> must coincide with the
  model / script filename stem. For example, a test case using ``case01.mod``
  only could be
  called ``case01 objno=5``.

* **"tags": ["linear", "continuous"]**: tags specifying model type, the case
  is executed only if the tags are a subset of the solver's ones. Except the
  tag **"run"** which means the test case is an AMPL script.

* *"files": ["diet.mod", "diet.dat"]*

* *"objective": value*: expected objective value.

* *"options": { "ANYSOLVER_options": "iisfind=1", "baron_options": "iisfind=12", "send_statuses": "0" }*.
  Option key ending with ``ANYSOLVER_options`` is for any solver, except when
  a solver-specific key is present (like ``baron_options``.)

* *"values": { "X[0].iis": "upp", ... }*. Expected values or expressions,
  in the form AMPL ``display`` command would accept.

  * For example, to check *logical expressions*, use if/then:

    .. code-block:: json

        "values": {
          "if color['Belgium'] != color['France'] then 1": 1,
          "if abs(x) < 1e-3 then 1": 1,
          "solve_result_num": 0
        }


.. _unit_tests:

Unit tests
----------

You can employ unit tests as well, see folder :file:`test`.
