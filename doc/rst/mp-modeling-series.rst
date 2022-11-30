
.. _mp-modeling-series:

MP Modeling Series: small examples
---------------------------------------

This series presents short discussions about common logical and
non-linear modeling tasks.

MP Modeling Series #1: a disjunction
****************************************

Assume the following constraint:
at most one of two variables *x*, *y* can be positive.
For the new MP Library-based drivers (gurobi, highs, copt),
as well as for Constraint Programming solvers (ilogcp, gecode, jacop),
this goes via AMPL logical operators:

.. code-block:: ampl

     x <= 0 or y <= 0.


Complete examples implementing this constraint in different ways:

1) With MP or CP, using the *or* operator:

   .. code-block:: ampl

        var x >= -1000 <= 1000;
        var y >= -1000 <= 1000;
        maximize obj: x + y;
        s.t. c1: x <= 0 or y <= 0;
        option solver gurobi;
        solve;
        display x, y;
        option solver highs;
        solve;
        display x, y;


2) With MP, using implication:

   .. code-block:: ampl

        var x >= -1000 <= 1000;
        var y >= -1000 <= 1000;
        maximize obj: x + y;
        s.t. c1: x > 0 ==> y <= 0;


3) Without MP:

   .. code-block:: ampl

        var x >= -1000 <= 1000;
        var y >= -1000 <= 1000;
        var xx binary;
        var yy binary;
        maximize obj: x + y;
        s.t. c1: x <= xx1000;
        s.t. c2: y <= yy1000;
        s.t. c3: xx+yy <= 1;
        option solver gurobiasl;
        solve;
        display x, y;
        option solver highs;
        solve;
        display x, y;

