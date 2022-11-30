
.. _numerical_accuracy:

Numerical accuracy
------------------------

Mathematical Programming solvers typically work with finite-precision numbers, which
leads to concerns on numerical stability.

Relational operators
******************************

The MP library simplifies relational operators into "indicator" constraints.
Solvers natively supporting indicators, usually handle them in a numerically stable way.
Otherwise, they have to be linearized by the so-called "big-M" constraints. The big-M
constants require finite bounds on expressions. For numerical stability these bounds should
not exceed the reciprocal of the integrality tolerance (option *inttol*). A default
big-M value can be set with the option *cvt:bigM*.

Piecewise-linear functions
*****************************

Piecewise-linear expressions can be modeled in AMPL directly, or arise from
approximations of other functions. Solvers which support PL expressions,
usually handle them algorithmically in a numerically stable way. Otherwise,
if PL expressions are linearized, it is recommended to have the argument
and result variables bounded in [-1e+4, 1e+4] (for approximated nonlinear functions,
hard bounds of up to [-1e+6, 1e+6] are imposed). The stability can be improved
in some cases by decreasing integer tolerance, Gurobi's *intfocus* and
*numfocus* options, switching off presolve in the solver, and other tuning measures.

