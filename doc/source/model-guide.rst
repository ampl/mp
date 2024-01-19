.. _modeling-guide:

=========================================
Modeling guide for MP-based AMPL solvers
=========================================


Ever wondered how to model logical and non-linear constraints? For example:

- At most one of the variables *x*, *y* can be positive

- All variables in a group should have different values

- Piecewise-linear approximation of *y = sin(x)*.

MP automates many of such modeling tasks by reformulating AMPL nonlinear
constraints suitably for the given solver.
A series of small modeling tasks like these are handled in the
`MP Modeling Series <https://discuss.ampl.com/c/modeling/42>`_,
while below follows a comprehensive guide.


.. toctree::
   :maxdepth: 2

   self
   modeling-overview
   modeling-expressions
   modeling-efficiency
   modeling-numerics
   modeling-tools
   modeling-troublesh
