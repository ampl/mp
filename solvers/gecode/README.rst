The solver ``gecode`` uses `Gecode constraint development environment
<http://www.gecode.org/>`_ to solve constraint programming problems.
Solver binaries are available for download from the
`Open Source Solvers <http://ampl.com/products/solvers/open-source#gecode>`_
page.

Normally the gecode solver is invoked by AMPL's ``solve`` command,
which gives the invocation
::

     gecode stub -AMPL

in which ``stub.nl`` is an AMPL generic output file (possibly written
by ``ampl -obstub`` or ``ampl -ogstub``).  After solving the problem,
the solver writes a ``stub.sol`` file for use by ampl's ``solve`` and
``solution`` commands. When you run ampl, this all happens automatically
if you give the AMPL commands
::

     option solver gecode;
     solve;

You can control the solver by setting the environment variable
``gecode_options`` appropriately (either by using ampl's ``option`` command,
or by using the shell's ``set`` and ``export`` commands before you invoke ampl).
You can put one or more (white-space separated) option assignments in
``$gecode_options``. The option ``version`` doesn't take a value:

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

     gecode -=

See `Gecode Options for AMPL
<http://ampl.com/products/solvers/gecode-options/>`_ for the full list of options.

solve_result_num values
=======================

Here is a table of ``solve_result_num`` values that ``ilogcp`` can return
to an AMPL session, along with the text that appears in the associated
``solve_message``.

=====   =================================================
Value   Message
=====   =================================================
    0   optimal solution (for an optimization problem) or
        feasible solution (for a satisfaction problem)
  200   infeasible problem
  400   time limit
  401   node limit
  402   fail limit
  403   solution limit
  600   interrupted
=====   =================================================

If you invoke ``gecode stub -AMPL`` or ``gecode stub``, you can also
supply additional command-line arguments of the form name=value.
Such arguments override specifications in ``$gecode_options``.  Example::

     ampl -obfoo foo.model foo.data
     nohup gecode -s foo 2>>err&

to solve a problem whose solution will take a while; after it finishes,
::

     ampl foo.model foo.data -
     solution foo.sol;
     display ... /* things involving the computed solution */;

(Here, ``-`` denotes standard input, and ampl reads the ``solution...``
and ``display...`` lines.)

Suffixes
========

You can use the suffix ``icl`` to specify integer consistency level for
constraints::

  subj to c1: alldiff ({i in 1..n} q[i]) suffix icl icl_dom;

This suffix and possible values for it are defined in ``gecode.ampl``.
Requires AMPL version 20130906 or later.

------------

See also `AMPL extensions for constraint programming
<http://ampl.com/resources/logic-and-constraint-programming-extensions/>`_.

If you have questions about or find bugs with this stuff,
please contact `Victor Zverovich <mailto:viz@ampl.com>`_.
