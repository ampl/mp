
.. _modeling-troublesh:

Troubleshooting
---------------------------

Always look for a more efficient model.
Bounds could be tightened; coefficients' magnitude can be reduced
by rescaling of the data; you might not need certain constraints
or logical conditions when using a different approach.
Sometimes, just a different
formulation can help, for example, a logical condition
can be manually linearized.
See sections :ref:`efficiency` and :ref:`numerical_accuracy`.

Play with `AMPL options <https://dev.ampl.com/ampl/options.html>`_
and :ref:`solver-options`.
Prominent ones are AMPL presolve
(switch off: ``ampl: option presolve 0;``) and solver's presolve
(``ampl: option gurobi_options 'presolve=0';``) and others
(tolerances, *numfocus*, *intfocus*, etc.)

To see what MP and/or the solver do with your model, export
the solver's received model, and, if possible, the solver's presolved model:

.. code-block:: ampl

    option gurobi_options 'writeprob=disj.lp writepresolved=disj_pre.lp';
    option gurobi_auxfiles rc;     ## To use var/con names
    solve;

MP offers ways to :ref:`explore automatic reformulations<explore-reformulations>`.

If you decide to contact AMPL or solver support, please provide a (possibly reduced)
version of your model reproducing the issue. Please also provide
the AMPL solver version by running ``highs -v``
or solver logs obtained with options ``version outlev=1``.
If you cannot reduce the model and don't
want to show it, give us just the NL file produced by

.. code-block:: ampl

    option auxfiles ''; write bdisj;

For solver support, use MPS or another full-precision model format.
