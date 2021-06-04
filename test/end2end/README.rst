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
  them in the local ``modellist.json`` (please refer to existing examples for
  the format.)
