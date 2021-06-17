AMPL/MP/END-2-END
=================

.. image:: https://travis-ci.org/ampl/mp.png?branch=master
  :target: https://travis-ci.org/ampl/mp

.. image:: https://ci.appveyor.com/api/projects/status/91jw051om9q8pwt9
  :target: https://ci.appveyor.com/project/vitaut/mp

An open-source solver testing library.


HOWTOs
------

* To **run all tests** for a given ``solver``, enter the command:

  .. code-block:: console
  
      python3 test/end2end/run.py solver
      
  The ``solver`` and ``ampl`` are expected to be on the system path.
  
* To run a **subset of the test cases**, ``cd`` into the corresponding
  subfolder of ``test/end2end/cases``, or use the ``--dir`` or
  ``--nonrecursive`` options.
  
* To **add a new solver** to the library, derive a corresponding class in
  ``Solver.py`` and list it in ``SolverCollection.py``.
  
* To **add new test cases**, add the model/data/AMPL script files and describe
  them in the local ``modellist.json`` having the following format. The top JSON
  object is an array of test cases. Each element is a dictionary with the
  following items, where non-compulsory items are italicized:
  
      * **"name": "<name>"**: <name> is the case name, must coincide with the
        model / script stem. For example, test case using ``diet.mod`` should be
        called ``diet``. There can be several cases with equal name.
      
      * **"tags": ["linear", "continuous"]**: tags specifying model type, the case
        is executed only if the tags are a subset of the solver's ones.
      
      * *"files": ["diet.mod", "diet.dat"]*
      
      * *"objective": value*: expected objective value.
      
      * *"options": { "ANYSOLVER_options": "iisfind=1", "baron_options": "iisfind=12", "send_statuses": "0" }*.
        Option key ending with ``SOLVER_options`` is for any solver, except when
        a solver-specific key is present (like ``baron_options``.)
        
      * *"values": { "X[0].iis": "upp", ... }*. Expected values.
