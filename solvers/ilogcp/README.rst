Ilogcp
======

The solver ``ilogcp`` connects IBM ILOG `CPLEX CP Optimizer
<http://www-01.ibm.com/software/integration/optimization/cplex-cp-optimizer/>`_
and `CPLEX <http://www-01.ibm.com/software/integration/optimization/cplex-optimizer/>`_
to AMPL via the Concert API.
It fully supports `AMPL extensions for logic and constraint programming
<http://ampl.com/resources/logic-and-constraint-programming-extensions/>`_ and
can handle wide range of problems including linear, mixed-integer, quadratic and
constraint programming problems.

Normally the ilogcp solver is invoked by AMPL's ``solve`` command, which
gives the invocation
::

     ilogcp stub -AMPL

in which ``stub.nl`` is an AMPL generic output file (possibly written
by ``ampl -obstub`` or ``ampl -ogstub``). After solving the problem,
the solver writes a ``stub.sol`` file for use by ampl's ``solve`` and
``solution`` commands. When you run ampl, this all happens automatically
if you give the AMPL commands
::

     option solver ilogcp;
     solve;

You can control the solver by setting the environment variable
``ilogcp_options`` appropriately (either by using ampl's ``option`` command,
or by using the shell's ``set`` and ``export`` commands before you invoke ampl).
You can put one or more (white-space separated) option assignments in
``$ilogcp_options``.  The option ``version`` doesn't take a value:

=======      ==================================================
Phrase       Meaning
=======      ==================================================
version      Report version details before solving the problem.
=======      ==================================================

Others are name-value pairs separated by '=', as in
::

     timelimit=600

which limits CP Optimizer search time to 600 seconds.  Options such
as ``logverbosity``, that take enumerated list of values, accept both numeric
and string values, so
::

     logverbosity=terse

is equivalent to
::

     logverbosity=1

In particular, switches that take values 0 or 1 also accept values
``off`` or ``on``.

The following command prints the full list of options with descriptions::

     ilogcp -=

See `IBM ILOG CPLEX CP Optimizer Options for AMPL
<http://ampl.com/products/solvers/ilogcp-options/>`_ for the full list of options.

solve_result_num values
-----------------------

Solve result codes follow MP conventions: 0 denotes an optimal solution
(or a solution to a satisfaction problem), 100 a feasible solution without
proven optimality, 200 infeasibility, and 300 unboundedness with a feasible
solution (350 without one).
Limits use 400-series codes, including time limits (402/472), solution
limits (408), interrupts (401/471), fail limits (430/490), branch limits
(431/491), and choice point limits (432/492). Infeasible or unbounded
is 450; errors use 500-series codes. Where paired codes are listed,
the first indicates a feasible solution and the second no feasible solution.

The CPLEX backend also reports its applicable limit and numerical statuses.
Run ``ilogcp -!`` for the full list of result codes and messages.

------------

If you invoke ``ilogcp stub -AMPL`` or ``ilogcp stub``, you can also
supply additional command-line arguments of the form name=value.
Such arguments override specifications in ``$ilogcp_options``.  Example::

     ampl -obfoo foo.model foo.data
     nohup ilogcp -s foo timing=1 2>>times&

to solve a problem whose solution will take a while; after it finishes,
::

     ampl foo.model foo.data -
     solution foo.sol;
     display ... /* things involving the computed solution */;

(Here, ``-`` denotes standard input, and ampl reads the ``solution...``
and ``display...`` lines.)
