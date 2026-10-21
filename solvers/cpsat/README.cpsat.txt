AMPL driver for OR-Tools CP-SAT
===============================

The solver "cpsat" connects Google OR-Tools CP-SAT to AMPL, using the MP library.

CP-SAT works with integer variables and integer coefficients only. Besides
linear constraints, the following are passed natively to the solver:
alldiff, count, numberof, if-then-else, min, max, abs, logical and/or,
and indicator constraints. Other constraints are reformulated by MP when possible.

To select the solver:

  ampl: option solver cpsat;

and to set options:

  ampl: option cpsat_options 'timelim=60 threads=4 outlev=1';

Use "cpsat -=" to list all options.
