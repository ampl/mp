.. _solution-check:


Solution check
---------------------

Solutions obtained from the solver are automatically checked
for correctness with given tolerances
(see driver options ``sol:chk:...``.)


"Realistic" solution check
******************************

In this mode, variable values are taken as they were reported by the solver
(with possible modifications using AMPL's options
``solution_round`` and ``solution_precision``, see driver options
``sol:chk:no...``.). This check is enough for most practical situations.

.. code-block:: ampl

		------------ WARNINGS ------------
		WARNING:  "Solution Check"
				 [ sol:chk:feastol=1e-06, sol:chk:inttol=1e-05,
				 solution_round='', solution_precision='' ]
		Algebraic expression violations:
			- 1 original expression(s) of type ':quadrange',
					up to 1E+00 (item 'socp[13]')


"Idealistic" solution check
******************************

In this mode, expressions are recomputed. Consider the following example.

.. code-block:: ampl

		var x >=0, <=100;
		maximize Total: if x<=5 and x>=5.00000000001 then 10;

Most solvers apply a constraint feasibility tolerance of the order :math:`10^{-6}`.

.. code-block:: ampl

		ampl: option solver gurobi;
		ampl: solve;
		Gurobi 10.0.2: optimal solution; objective 10
		0 simplex iterations

		------------ WARNINGS ------------
		WARNING:  "Solution Check (Idealistic)"
				 [ sol:chk:feastol=1e-06, sol:chk:inttol=1e-05,
				 solution_round='', solution_precision='' ]
		Objective value violations:
			- 1 objective value(s) violated,
					up to 1E+01

		ampl: display x;
		x = 5

We see that ``x=5`` satisfies the ``if`` with that tolerance.
Thus, our realistic check passes, but the idealistic check complains.
Indeed, if we ask AMPL to recompute the objective value:

.. code-block:: ampl

		ampl: display Total;
		Total = 0

we see that AMPL does it "idealistically"
(it does not know about solver tolerances.)

To see which expressions cause the violation,
use driver option ``chk:mode``:

.. code-block:: ampl

		ampl: option gurobi_options 'chk:mode=1023';
		ampl: solve;
		Gurobi 10.0.2:   sol:chk:mode = 1023
		Gurobi 10.0.2: optimal solution; objective 10
		0 simplex iterations

		------------ WARNINGS ------------
		WARNING:  "Solution Check (Idealistic)"
				 [ sol:chk:feastol=1e-06, sol:chk:inttol=1e-05,
				 solution_round='', solution_precision='' ]
		Algebraic expression violations:
			- 1 original expression(s) of type ':ifthen',
					up to 1E+01
		Logical expression violations:
			- 1 original expression(s) of type ':and'
		Objective value violations:
			- 1 objective value(s) violated,
					up to 1E+01

*Hint*: to display AMPL model names,
set ``option (solver_)auxfiles rc;`` as follows:

.. code-block:: ampl

		ampl: option gurobi_auxfiles rc;
		ampl: solve;
		Gurobi 10.0.2:   sol:chk:mode = 1023
		Gurobi 10.0.2: optimal solution; objective 10
		0 simplex iterations

		------------ WARNINGS ------------
		WARNING:  "Solution Check (Idealistic)"
				 [ sol:chk:feastol=1e-06, sol:chk:inttol=1e-05,
				 solution_round='', solution_precision='' ]
		Algebraic expression violations:
			- 1 original expression(s) of type ':ifthen',
					up to 1E+01 (item 'Total_11_')
		Logical expression violations:
			- 1 original expression(s) of type ':and'
					(item 'Total_7_')
		Objective value violations:
			- 1 objective value(s) violated,
					up to 1E+01 (item 'Total')
