
.. _howto-test:

Testing
=======

Choose the check that exercises the behavior you changed.

.. list-table::
   :header-rows: 1
   :widths: 17 33 29 21

   * - Check
     - What it exercises
     - What it needs
     - How it is run
   * - :ref:`CTest <unit_tests>`
     - MP library code, converters and NL Writer behavior.
     - A build with ``BUILD_TESTS=ON``. Optional tests may need their
       component or solver SDK.
     - Locally or in jobs that build and run the registered CTest targets.
   * - :ref:`CLI end-to-end <run_end2end_tests>`
     - The model, driver and solver together, including status, objective,
       suffixes and residuals when specified by the case.
     - A driver executable, its runtime and license, Python dependencies
       and AMPL for cases written as AMPL models.
     - Run ``test/end2end/run.py`` locally or in a solver test job.
   * - :ref:`Windows AMPLS <ampls_repeat_test>`
     - Reporting and solving again on the same native model.
     - Windows Python, AMPL, a shared driver and its matching solver DLLs.
     - Manual. The default CI workflow does not invoke it.

Check the selected CI job to see which suites it invokes and whether
their failures fail the job. Publishing results alone does not establish
that test failures are enforced. The standalone MP wheel-building
workflow does not run the driver suites above.

:ref:`Add a regression case <add_new_test_cases>` when changing an
essential driver behavior.


End-to-end tests
----------------

End-to-end tests have proven most important.


.. _run_end2end_tests:

Run end-to-end tests
~~~~~~~~~~~~~~~~~~~~

* To **run all tests** for one or several ``solver``\ s, enter the command:

  .. code-block:: console

      python3 test/end2end/run.py solver [another_solver [...]] [--options 'acc:max=0']

  The ``solver`` and ``ampl`` executables are expected to be on the system path,
  or use options ``--ampl`` (default ``$MP_TEST_AMPL``), ``--bin_path``.
  With ``--ampl_fallback`` (default ``$MP_TEST_AMPL_FALLBACK``), a release
  AMPL is used instead when the selected one is not found or is a sanitizer
  build, whose diagnostics would appear in every solver log.

* Detailed results are saved into a CSV report file, see ``--reportstub``.

* To run a **subset of the test cases**, ``cd`` into the corresponding
  subfolder of ``test/end2end/cases``, or use the ``--dir`` or
  ``--nonrecursive`` options.

  For example, run the ``fast`` folder with one solver thread and a
  five-second per-instance limit::

      python3 test/end2end/run.py highs --bin_path path/to/bin --dir test/end2end/cases/categorized/fast --timeout 5 --nthreads 1 --reportstub fast-highs

  Check the final summary for solver failures, aborted runs, AMPL/Python
  failures and skips. The runner can exit successfully and still report
  failed comparisons. Skipped cases did not exercise the driver.

  For a regression comparison, run the baseline and candidate with the
  same AMPL and solver SDK versions, case set, options, thread count and
  timeout. Use distinct report stubs and compare by folder and case name.
  Record added and cleared failures separately from timeouts and skips.
  If a case unexpectedly skips, check its feature tags and solver
  restrictions before drawing a coverage conclusion.


.. _ampls_repeat_test:

Repeated solves through AMPLS on Windows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ordinary end-to-end runner starts a fresh driver invocation for each
solve. The Windows integration test
:file:`test/solvers/ampls-fixed-model-test.py` exercises a different path:
it loads a model once, solves and reports it with ``mip:basis=1``, reverses
the objective sense on the same native problem, then solves and reports
again. It checks both fixed solves, the unchanged native problem, and the
change in the integer solution. It uses the CPLEX MIQCP case in
``fast/mip_basis`` and the Xpress ``nl_dual_01`` case in ``fast/suf_std``.

Run this test manually from the MP source root. It requires Windows
Python, a licensed AMPL, and the selected shared driver and matching
solver SDK DLLs. The examples use CPLEX 22.2 and Xpress::

    python test/solvers/ampls-fixed-model-test.py cplex C:/path/to/bin/cplex-lib.dll C:/path/to/bin/cplex2220.dll --ampl C:/path/to/ampl.exe
    python test/solvers/ampls-fixed-model-test.py xpress C:/path/to/bin/xpress-lib.dll C:/path/to/bin/xprs.dll --ampl C:/path/to/ampl.exe

Use the matching DLL filename for your SDK version and quote paths with
spaces. Other runtime DLLs must be available to the loader. This test runs
outside CTest and the standard end-to-end runner. The default CI
configuration does not run it.

Compulsory configurations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following test configurations are compulsory for a release.

* All important/representative solvers. While some tests use `acc:...` options
  to enforce linearization etc., open-source solvers use it everywhere.
  If desired, such options can be passed to all test cases using the `--options` switch,
  e.g., `--options 'acc:indle=0 acc:indge=0 acc:indeq=0 acc:and=0 acc:or=0'`.

* Using solution checker with `--options "chk:fail chk:feastol=1e-2[ chk:mode=1023]"`.


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

To **add new test cases**, add the model/data/AMPL script/NL files in
a subfolder of :file:`test/end2end/cases/categorized/` and describe
them in the local ``modellist.json`` having the following format.


Test specification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In ``modellist.json``, the top JSON
object is an array of test cases. Each element is a dictionary with the
following items, where non-compulsory items are italicized:

* **"name": "<name>"**: case name. Unless *"files"* is present
  (see below), the first word of <name> must coincide with the
  model / script filename stem. For example, a test case using ``case01.mod``
  only could be
  called ``case01 objno=5``.

* **"tags": ["linear", "continuous"]**: tags specifying model type, the case
  is executed only if the tags are a subset of the solver's ones. Except the
  tag **"script"** which means the test case is an AMPL script. Such a
  case executes its own solves and runs in the file's directory, so relative
  model, data, and include paths are resolved there.
  ``--timeout`` also protects script execution, including its solves.
  List **"continuous"** whenever the model has continuous variables
  (including AMPL's ``var x in 1..9``): drivers for integer models only
  skip such cases. Tag **"mpdriver"** marks a case that needs an MP driver
  (MP options or MP's handling of invalid options).

* *"files": ["diet.mod", "diet.dat"]*

* *"objective": value*: expected objective value.

* *"options": { "ANYSOLVER_options": "iisfind=1", "baron_options": "iisfind=12", "send_statuses": "0" }*.
  Option key ending with ``ANYSOLVER_options`` is for any solver, except when
  a solver-specific key is present (like ``baron_options``).
  Other keys starting with ``ANYSOLVER_`` set the solver's AMPL option
  of that suffix, e.g., ``"ANYSOLVER_auxfiles": "rc"`` sets
  ``highs_auxfiles`` for HiGHS (pass variable and constraint names).
  Moreover, ``mp_options`` applies to all MP drivers.
  ASL drivers skip cases with ``mp_options`` or with MP (namespaced)
  options such as ``cvt:bigm=1e5`` in their options; an ASL driver built
  with the MP2NL wrapper runs them through MP2NL (solver classes
  ``baron``, ``conopt``, ``knitro``, ``lindoglobal``, ``gurobiaslmp2nl``,
  ``xpressaslmp2nl``, ``cplexaslmp2nl``): ``mp2nl=1`` is passed in the
  driver's options, and the case's options go to ``mp2nl_options``.

* *"values": { "X[0].iis": "upp", ... }*. Expected values or expressions,
  in the form AMPL ``display`` command would accept. Only available for AMPL
  models / scripts. Example:

  .. code-block:: json

      "values": {
        "if color['Belgium'] != color['France'] then 1": 1,
        "if (forall {s in NS} (VAL_U_F2[s] = 1  ==>  U_F2[s] >= 0.0001  else  U_F2[s] = 0)) then 1": 1,
        "abs(x)": { "max": 1e-3 },
        "solve_result_num": 0
      }

  In particular:

  * The *expected value* can be a dictionary of the form `{ "min": ..., "max": ... }`:

    .. code-block:: json

        "max{i in 1 .. _nccons} abs(_ccon[i])": { "max": 1e-5 },
        "min{i in 1 .. _ncons} _con[i].slack": { "min": -1e-5 }

    This checks all complementarity and algebraic constraints.

  * To check *logical expressions*, use if/then:

    .. code-block:: json

        "if forall {i in 1.._nlogcons} _logcon[i] then 1": 1,

    This checks all logical constraints.


* *"output": ["Presolved model has 500 variables", "RHS range:  [1e-2, 1e4]"]*.
  Output chunks expected in the solver / AMPL log.


Checking reformulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to check characteristics of the models received by an MP
driver from the NL file, as well as of the model submitted to the solver.
For that, add option ``[tech:]debug=1``. Then, the driver outputs
problem suffixes of the form
``flat<N>__[<C>__]<itemtype>`` giving the number of model items
of various types and statuses. The specifiers:

* ``<N>`` is 0 for the initial flat model obtained from the NL file,
  1 for the reformulated model.

* If present, ``<C>`` can be ``A`` for the total number of items,
  ``U`` for the number of total items still relevant, but possibly reformulated.
  If not present, it is the number of active items at the end of step ``<N>``.

* ``<itemtype>`` can be a constraint type name, variable type, or objective type.

Run an AMPL session to see all such suffixes, or check existing test cases.

These counts depend on which constraints a driver accepts natively,
i.e., does not reformulate. Test drivers are therefore assigned a
*solver category* (``SolverCategory`` in the test scripts):

* ``MILP``: linear constraints only (HiGHS, GCG, BARON);
* ``MIPWithIndicators``: plus indicator constraints (CPLEX, COPT, MOSEK);
* ``MILPWithLogicals``: plus ``and``, ``or`` (Gurobi, Xpress, SCIP);
* ``CP``: plus ``count``, ``numberof``, ``alldiff``, if-then (CP-SAT);
* ``CPLegacy``: drivers on the old MP framework, without the model converter
  and its counts (Gecode, ilogcp);
* ``MP2NL``: the MP2NL driver, writing the NL format; for example,
  complementarity is passed as ``nlcompl`` (unless ``acc:nlcompl=0``).
  Set ``debug=1`` in *mp2nl_options*: with *mp_options*, an MP subsolver
  reports its own counts, overriding those of MP2NL.

Nonlinear functions do not influence the category, except for ``MP2NL``.
The test scripts warn when the ``acc:`` options of a driver imply a category
other than the declared one.
A case gives *"values"* as for one category and corrections for others in
*"values_per_solver_category"*; ``null`` drops a check:

.. code-block:: json

    "values": {
      "Initial.flat1__U__indeq": 6,
      "Initial.flat1__linle": 14
    },
    "values_per_solver_category": {
      "CP": { "Initial.flat1__U__indeq": 3, "Initial.flat1__linle": 8 },
      "CPLegacy": { "Initial.flat1__U__indeq": null, "Initial.flat1__linle": null }
    }

Drivers without a category use *"values"* only.


.. _unit_tests:

Unit tests
----------

Unit tests are in the folder :file:`test`.
Compulsory are tests `converter-flat-test` and `converter-mip-test`.
Configure MP with ``BUILD_TESTS=ON``. Use ``NLW2_LIB=ON`` to include the
NL Writer library and its solver-call tests when examples are disabled.
For a multi-configuration build, run ``ctest --test-dir path/to/build -C
Release --output-on-failure`` after building the test targets.

With the end-to-end Python dependencies installed, run the script-timeout
regression checks with ``python3 -m unittest discover -s
test/end2end/scripts/python -p test_script_timeout.py``. They include a
short-lived subprocess test of solver-descendant cleanup.
