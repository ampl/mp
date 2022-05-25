Reusable components and solver driver API
=========================================


This section highlights the major library components,
ranging from basic building blocks necessary for any AMPL driver,
to more advanced utilities.


NL and SOL files
----------------

Any AMPL solver driver currenlty needs to input
an NL file and report results in a SOL file.
Corresponding APIs are described below.

NL file reader
~~~~~~~~~~~~~~

MP provides a high-performance :doc:`.nl file reader <nl-reader>`
which is up to `6x faster
<http://zverovich.net/slides/2015-01-11-ics/socp-reformulation.html#/14>`_
than the one provided by the traditional
`AMPL Solver Library (ASL)
<https://ampl.com/resources/learn-more/hooking-your-solver-to-ampl/>`_.


NL format
^^^^^^^^^

`NL <https://en.wikipedia.org/wiki/Nl_(format)>`_ is a format for representing
optimization problems in discrete or continuous variables. It is described in
the technical report `Writing .nl Files <https://ampl.github.io/nlwrite.pdf>`_.

The NL format supports a wide range of problem types including but not limited
to the following areas of optimization:

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

This section describes the C++ API of an NL reader which is

* Reusable: the reader can be used to process NL files in different ways
  and not limited to a single problem representation
* High performance: fast `mmap <http://en.wikipedia.org/wiki/Mmap>`_-based reader
  with `SAX <http://en.wikipedia.org/wiki/Simple_API_for_XML>`_-like API and no
  dynamic memory allocations in the common case
* Easy to use: clean, modern code base and simple API
* Complete: supports all NL constructs including extensions implemented in
  AMPL Solver Library
* Reliable: extensively and continuously tested on a variety of platforms


Easy-to-use functions
^^^^^^^^^^^^^^^^^^^^^

The ``mp/nl.h`` header only contains declarations of
`mp::ReadNLFile` and `mp::ReadNLString`, and can be used to read the standard optimization problem
object of class `mp::Problem`, for example:

.. code-block:: c++

   #include "mp/nl.h"
   #include "mp/problem.h"

   mp::Problem p;
   ReadNLFile("diet.nl", p);


Full NL-reader API
^^^^^^^^^^^^^^^^^^

If you want to provide a custom NL handler, include ``mp/nl-reader.h`` instead.
Class `mp::NLHandler` can be customized for most efficient translation of NL format into
solver API using the `mp::ProblemBuilder` concept.

* Alternatively, convenience classes `mp::Problem` and `mp::ColProblem`
  can be used for intermediate storage of the NL model. From there,
  `mp::ExprVisitor` and `mp::ProblemFlattener` walk the NL forest.

A few examples are in
`nl-example.cc <https://github.com/ampl/mp/blob/master/src/nl-example.cc>`_, `mp::Problem`,
`SCIP NL file reader <https://scipopt.org/>`_.
Note that ``mp/nl.h`` is a much smaller header than ``mp/nl-reader.h`` so prefer
it unless you need access to the full NL reader API, described below.


* `mp::ReadNLFile`, `mp::ReadNLString` read NL model

* `mp::NLHandler`, `mp::NullNLHandler` provide interface for a custom NL handler

* `mp::NLHeader` stores problem information

* `mp::ReadError`, `mp::BinaryReadError`

* `mp::arith::Kind`

* `mp::READ_BOUNDS_FIRST` can be passed as a flag to `mp::ReadNLFile`

* `mp::MAX_AMPL_OPTIONS` is the maximum number of options reserved for AMPL use in NL and SOL formats



SOL file writer
~~~~~~~~~~~~~~~

A minimal .sol file writer can be implemented by parameterizing the
`internal::AppSolutionHandlerImpl` (or `internal::SolutionWriterImpl`)
templates by minimal implementations of the `mp::BasicSolver` and
`mp::ProblemBuilder` interfaces.


Driver logic
------------

Using the `mp::Backend` and the derived classes is now the
recommended approach to building a new solver interface.
They provides a convenient API for common solver driver actions,
options and suffixes.


Backend, MIPBackend
~~~~~~~~~~~~~~~~~~~

* `mp::Backend`, `mp::MIPBackend` standardize some common AMPL app behaviour, such as
  solver messages and status reporting, simplex basis statuses, and suffix I/O


Solver, SolverImpl [deprecated]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* `mp::Solver` and `mp::SolverImpl` enable very basic standard behaviour
  (e.g., multiobj, solution output)


* Classes `mp::BasicSolver`, `mp::Backend` and `mp::MIPBackend`
  standardize solver behavior such as common options and suffixes
  and are recommended for new interfaces.



Model/solution conversions
--------------------------

* Class `mp::BasicModelManager` standardizes the workflow of
  model input and results output.

* Classes `mp::FlatConverter` and `mp::MIPFlatConverter` facilitate conversion of
  NL expressions which are not natively accepted by a solver into simpler forms.

  * :ref:`Logical and CP constraints <modeling-guide>` are supported.

* Class `mp::pre::Presolver` pre- and postsolves solutions and suffixes.



C++ ASL adapter
---------------

An efficient type-safe `C++ adapter for the traditional ASL library
<https://github.com/ampl/mp/tree/master/src/asl>`_ for
connecting solvers to AMPL and other systems. ASL has many additional functions,
such as writing NL files and automatic differentiation.




More details
------------

This section overviews the common API components in more detail.

For a complete API reference, see the :ref:`index <genindex>`.


Problem builders
~~~~~~~~~~~~~~~~

* `mp::ProblemBuilder`

* `mp::ColProblem`


Problem representation
~~~~~~~~~~~~~~~~~~~~~~

A standard representation of a model, convenient for intermediate storage.
Can be converted into solver API by a subclassed `mp::ExprVisitor`.

* `mp::ProblemInfo`, `mp::var::Type`, `mp::obj::Type`, `mp::func::Type`, `mp::ComplInfo`

* `mp::Problem`, `mp::BasicProblem`


Expression forest walkers
~~~~~~~~~~~~~~~~~~~~~~~~~

Typesafe expression walkers for models stored in memory.

* `mp::expr::Kind`, `mp::expr::str`, `mp::expr::nl_opcode`

* `mp::BasicExprVisitor`, `mp::ExprVisitor`, `mp::ExprConverter`

* `mp::ProblemFlattener`


Model and solution management
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Class `mp::BasicModelManager` standardizes the workflow of
  model input and results output.


Model conversion and presolve
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Classes `mp::FlatConverter` and `mp::MIPFlatConverter` facilitate conversion of
  NL expressions which are not natively accepted by a solver into simpler forms.

  * `mp::ConstraintKeeper` stores constraints in FlatConverter.

* Class `mp::pre::Presolver` pre- and postsolves solutions and suffixes.



Standard AMPL driver logic
~~~~~~~~~~~~~~~~~~~~~~~~~~

* `mp::Backend`, `mp::MIPBackend`

* `mp::Solver`, `mp::SolverImpl` [deprecated]


Solution status
~~~~~~~~~~~~~~~

* `mp::sol::Status`


Suffixes
~~~~~~~~

* `mp::suf::Kind`, `mp::SuffixDef`

* Standard suffix value enums: `mp::IISStatus`, `mp::BasicStatus`


