Jacop
=====

The solver ``jacop`` uses `Java Constraint Programming solver (JaCoP)
<http://jacop.osolpro.com/>`_ to solve constraint programming problems.
It fully supports `AMPL extensions for constraint programming
<http://www.ampl.com/NEW/LOGIC>`_. Solver binaries are available for
download from the `Open Source Solvers
<http://ampl.com/products/solvers/open-source#jacop>`_ page.

Normally the jacop solver is invoked by AMPL's ``solve`` command,
which gives the invocation
::

     jacop stub -AMPL

in which ``stub.nl`` is an AMPL generic output file (possibly written
by ``ampl -obstub`` or ``ampl -ogstub``).  After solving the problem,
the solver writes a ``stub.sol`` file for use by ampl's ``solve`` and
``solution`` commands. When you run ampl, this all happens automatically
if you give the AMPL commands
::

     option solver jacop;
     solve;

You can control the solver by setting the environment variable
``jacop_options`` appropriately (either by using ampl's ``option`` command,
or by using the shell's ``set`` and ``export`` commands before you invoke ampl).
You can put one or more (white-space separated) option assignments in
``$jacop_options``. The option ``version`` doesn't take a value:

=======      ==================================================
Phrase       Meaning
=======      ==================================================
version      Report version details before solving the problem.
=======      ==================================================

Others are name-value pairs separated by '=', as in
::

     timelimit=600

which limits solution time to 600 seconds.

The following command prints the full list of options with descriptions::

     jacop -=

solve_result_num values
-----------------------

Here is a table of ``solve_result_num`` values that ``jacop`` can return
to an AMPL session, along with the text that appears in the associated
``solve_message``.

=====   =================================================
Value   Message
=====   =================================================
    0   optimal solution (for an optimization problem) or
        feasible solution (for a satisfaction problem)
  100   feasible solution
  200   infeasible problem
  400   limit
  600   interrupted
=====   =================================================

------------

If you invoke ``jacop stub -AMPL`` or ``jacop stub``, you can also
supply additional command-line arguments of the form name=value.
Such arguments override specifications in ``$jacop_options``.  Example::

     ampl -obfoo foo.model foo.data
     nohup jacop -s foo 2>>err&

to solve a problem whose solution will take a while; after it finishes,
::

     ampl foo.model foo.data -
     solution foo.sol;
     display ... /* things involving the computed solution */;

(Here, ``-`` denotes standard input, and ampl reads the ``solution...``
and ``display...`` lines.)
