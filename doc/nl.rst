.. cpp:namespace:: mp

Nl Reader
=========

Nl is a format for representing optimization problems in discrete or continuous
variables. It is described in the technical report `Writing .nl Files
<https://cfwebprod.sandia.gov/cfdocs/CompResearch/docs/nlwrite20051130.pdf>`_.

The `nl format <https://en.wikipedia.org/wiki/Nl_(format)>`_ supports a wide
range of problem types including but not limited to the following areas of
optimization:

* `Linear programming <http://en.wikipedia.org/wiki/Linear_programming>`_
* `Quadratic programming <http://en.wikipedia.org/wiki/Quadratic_programming>`_
* `Nonlinear programming <http://en.wikipedia.org/wiki/Nonlinear_programming>`_
* `Mixed-integer programming <http://en.wikipedia.org/wiki/Linear_programming#Integer_unknowns>`_
* Mixed-integer quadratic programming with or without convex quadratic constraints
* Mixed-integer nonlinear programming
* `Second-order cone programming <http://en.wikipedia.org/wiki/Second-order_cone_programming>`_
* `Global optimization <http://en.wikipedia.org/wiki/Global_optimization>`_
* `Semidefinite programming <http://en.wikipedia.org/wiki/Semidefinite_programming>`_
  problems with bilinear matrix inequalities
* `Complementarity problems <http://en.wikipedia.org/wiki/Complementarity_theory>`_
  (MPECs) in discrete or continuous variables
* `Constraint programming <http://en.wikipedia.org/wiki/Constraint_programming>`_

This section describes the C++ API of an nl reader which is

* Reusable: the reader can be used to process nl files in different ways
  not limited to a single problem representation
* High performance: fast `mmap <http://en.wikipedia.org/wiki/Mmap>`_-based reader
  with `SAX <http://en.wikipedia.org/wiki/Simple_API_for_XML>`_-like API and no
  dynamic memory allocations in the common case
* Easy to use: clean, modern code base and simple API
* Complete: supports all nl constructs including extensions implemented in
  AMPL Solver Library
* Reliable: extensively and continuously tested on a variety of platforms

.. doxygenfunction:: ReadNLString(fmt::StringRef, Handler &, fmt::StringRef)

.. doxygenfunction:: ReadNLFile(fmt::StringRef, Handler &)
