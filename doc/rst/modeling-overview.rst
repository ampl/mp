.. _modeling-overview:

Modeling overview
-----------------


AMPL's newly extended C++ solver interface library, MP, is publicly
available in the `ampl/mp <https://github.com/ampl/mp>`_ repository.
Solver interfaces built with MP are able to handle a significantly
expanded range of model expressions.
Currently available MP-based solvers include:

- `gurobi <https://github.com/ampl/mp/tree/develop/solvers/gurobi>`_,
  an enhanced interface to the `Gurobi solver <https://ampl.com/products/solvers/solvers-we-sell/gurobi/>`_

- `copt <https://github.com/ampl/mp/tree/develop/solvers/copt>`_,
  an interface to `Cardinal Optimizer <https://ampl.com/products/solvers/solvers-we-sell/copt/>`_

- `highs <https://github.com/ampl/mp/tree/develop/solvers/highsmp>`_,
  an interface to the open-source `HiGHS solver <https://highs.dev/>`_ solver

- `cbc <https://github.com/ampl/mp/tree/develop/solvers/cbcmp>`_,
  an enhanced interface to the `CBC solver <https://ampl.com/products/solvers/open-source-solvers/>`_

- `xpress <https://github.com/ampl/mp/tree/develop/solvers/xpress>`_,
  an interface to `FICO Xpress <https://ampl.com/products/solvers/solvers-we-sell/xpress/>`_

Binaries for these solvers can be downloaded, in distribution
bundles and individually, through the `AMPL Portal <https://portal.ampl.com>`_.
More solvers will be added.



The expanded MP solver interface library offers new support
for the following categories of operators and expressions:

- Conditional operators: ``if-then-else``; ``==>``, ``<==``, ``<==>``
- Logical operators: ``or``, ``and``, ``not``; ``exists``, ``forall``
- Piecewise linear functions: ``abs``; ``min``, ``max``; ``<<breakpoints; slopes>>``
- Counting operators: ``count``; ``atmost``, ``atleast``, ``exactly``; ``numberof``
- Relational and comparison operators: ``>(=)``, ``<(=)``, ``(!)=``; ``alldiff``
- Complementarity operator: ``complements``
- Nonlinear operators and functions: ``*``, ``/``, ``^``; ``exp``, ``log``;
  ``sin``, ``cos``, ``tan``; ``sinh``, ``cosh``, ``tanh``
- Set membership operator: ``in``

Details and examples are given in the *Expressions supported* section below.
See also the individual solvers' documentation for details of solver-specific features:

- Choice between linearization in the interface and native solver support for some operations
- Handling of AMPL suffixes on constraints that are transformed by the interface

The slides from our presentation on
`Advances in Model-Based Optimization <https://ampl.com/MEETINGS/TALKS/2022_07_Bethlehem_Fourer.pdf>`_
provide overview of the MP interface library in the context of AMPL applications,
including comments on implementation and efficiency issues.

