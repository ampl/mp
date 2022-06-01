Reusable components and solver driver API
=========================================


This section highlights the major library components,
ranging from basic building blocks necessary for any AMPL driver,
to more advanced utilities.


.. _NL-SOL-files:

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


NL format
^^^^^^^^^

`NL <https://en.wikipedia.org/wiki/Nl_(format)>`_ is a format for representing
optimization problems in discrete or continuous variables. The format provides
linear constraints, as well as non-linear expression trees. It is described in
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
solver API using the :ref:`ProblemBuilder concept <problem-builders>`.

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

Writing solution/results output is easiest as part of the general workflow,
see :ref:`model-manager`.

A standalone .sol file writer could be implemented by parameterizing the
`mp::internal::AppSolutionHandlerImpl` (or `mp::internal::SolutionWriterImpl`)
templates by minimal implementations of the `mp::BasicSolver` and
`mp::ProblemBuilder` interfaces.


Driver logic
------------

Using the :ref:`mp::Backend and the derived classes <backend-classes>` is now the
recommended approach to building a new solver interface.
They provide a convenient API for common solver driver actions,
options and suffixes.
The high-level application structure is suggested as follows:

- :ref:`backend-app` --> :ref:`Custom Backend <backend-classes>` --> Solver.

Creating such driver from a template is
:ref:`described in the HowTo <howto-create-new-driver>`.


.. _backend-app:

BackendApp
~~~~~~~~~~

`mp::BackendApp` supports basic application functions, such as screen output
and interrupts handling. It calls a CustomBackend which should implement
the `mp::BasicBackend` interface.


.. _backend-classes:

The Backend classes
~~~~~~~~~~~~~~~~~~~

`mp::Backend` and `mp::MIPBackend` implement the `mp::BasicBackend` interface and
standardize some common AMPL app behaviour, such as
solver messages and status reporting,
simplex basis statuses, and other suffix I/O.
Their solver-specific subclasses can be customized for a particular solver.
They rely on the :ref:`model-manager` interface
for model and solution I/O.

As an example, if the driver should read and write simplex basis status suffixes,
the derived Backend class can declare

.. code-block:: c++

    ALLOW_STD_FEATURE( BASIS, true )
    SolutionBasis GetBasis() override;
    void SetBasis(SolutionBasis ) override;

and define the `GetBasis`, `SetBasis` methods.
See the classes' documentation
for further details.


.. _solver-classes:

Solver, SolverImpl [deprecated]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Classes `mp::SolverApp`, `mp::Solver` and `mp::SolverImpl` enable very basic
standard behaviour (e.g., multiobj, solution output). They are deprecated
in favor of the :ref:`BackendApp/Backend classes <backend-classes>` and
can be discontinued in future.



Model/solution I/O and conversions
----------------------------------

The tools presented in  this section standardize
model/solution I/O
(currently relying on :ref:`NL file input and SOL file output <NL-SOL-files>`)
and conversion for a particular solver.

.. _model-manager:

Model Manager
~~~~~~~~~~~~~

Class `mp::BasicModelManager` standardizes the interface for
model input and results output. This interface is used by the
:ref:`Backend classes <backend-classes>`.

* Current suggested implementations rely on `mp::ModelManagerWithProblemBuilder`
  which uses :ref:`NL file input and SOL file output <NL-SOL-files>` as well as
  a model converter. The model converter should implement the `mp::BasicConverter`
  interface and provide a :ref:`Problem Builder <problem-builders>`.


.. _problem-builders:

Problem builders
~~~~~~~~~~~~~~~~

Basic :ref:`Model/solution I/O <NL-SOL-files>` and
:ref:`model managers <model-manager>` rely on a `mp::ProblemBuilder` concept.

* A custom builder can pass the NL model directly into the solver. A few examples are in
  `nl-example.cc <https://github.com/ampl/mp/blob/master/src/nl-example.cc>`_, `mp::Problem`,
  `SCIP 8.0 NL file reader <https://scipopt.org/>`_.

* Alternatively, standard classes `mp::Problem` and `mp::ColProblem` provide intermediate
  storage for a problem instance. From `mp::Problem`,
  :ref:`conversion tools <problem-converters>`
  can be customized to transform the instance for a particular solver.


.. _problem-converters:

Problem converters
~~~~~~~~~~~~~~~~~~

Given a problem instance in the standard format `mp::Problem`, several
tools can be adapted to convert the instance for a particular solver.

* For :ref:`'flat' (expression-less) solvers <flat-solvers>`,
  `mp::ProblemFlattener` can walk the NL forest, passing flattened expressions as
  constraints to :ref:`flat-converters`. In turn, these
  facilitate conversion of flat constraints which are not natively accepted by a
  solver into simpler forms.

* For :ref:`expression-tree supporting solvers <expression-solvers>`,
  `mp::ExprVisitor` and `mp::ExprConverter` are efficient type-safe templates
  which can be customized to transform instances for a particular expression-based
  solver API.


.. _flat-converters:

Flat Model Converters
~~~~~~~~~~~~~~~~~~~~~

`mp::FlatConverter` and `mp::MIPFlatConverter`
facilitate conversion of flat models (i.e., models without expression trees).
Constraints which are not natively accepted by a
solver, are transformed into simpler forms. `mp::FlatConverter` and its subclasses
can be flexibly parameterized for a particular solver, preferably
via the solver's modeling API wrapper:

* `mp::BasicFlatModelAPI` is the interface via which `mp::FlatConverter` addresses
  the underlying solver. A subclassed wrapper, such as `mp::GurobiModelAPI`,
  signals its accepted constraints and which model conversions are preferable.
  For example, `GurobiModelAPI` declares the following in order to natively
  receive the logical OR constraint:

  .. code-block:: c++

      ACCEPT_CONSTRAINT(OrConstraint,
        Recommended,                       // Driver recommends native form
        CG_General)                        // Constraint group for suffix exchange
      void AddConstraint(const OrConstraint& );

  By default, if the driver does not mark a constraint as acceptable,
  `mp::FlatConverter` and its descendants attempt to simplify it. See the
  classes' documentation for further details.

* :ref:`value-presolver` transforms solutions and suffixes between the
  original NL model and the flat model.


.. _value-presolver:

Value Presolver
~~~~~~~~~~~~~~~

Class `mp::pre::ValuePresolver` manages transformations of solutions and suffixes
between the original NL model and the converted model. For driver architectures
with :ref:`model-manager`, the value presolver object must be shared between
the model converter and the :ref:`Backend <backend-classes>` to enable
solution/suffix transformations corresponding to those on the model, see
`mp::CreateGurobiModelMgr` as an example.


Invocation API
^^^^^^^^^^^^^^

To use the ValuePresolver API, the following classes are needed:

- `mp::pre::BasicValuePresolver` defines an interface for `mp::pre::ValuePresolver`.

- `mp::pre::ValueNode` temporarily stores values corresponding to a single type of
  model item (variables, constraints, objectives).

- `mp::pre::ValueMap` is a map of node values, where the key usually corresponds to
  an item subcategory. For example, Gurobi distinguishes attributes for the
  following constraint categories: linear, quadratic, SOS, general. Thus, the
  conversion graph needs to have these four types of target nodes for constraint
  values:

  .. code-block:: c++

    pre::ValueMapInt GurobiBackend::ConsIIS() {
      ......
      return {{{ CG_Linear, iis_lincon },
        { CG_Quadratic, iis_qc },
        { CG_SOS, iis_soscon },
        { CG_General, iis_gencon }}};
    }


- `mp::pre::ModelValues` is a class joining value maps for variables, constraints,
  and objectives. It is useful when the conversions connect items of different types:
  for example, when converting an algebraic range constraint to an equality
  constraint with a bounded slack variable, the constraint's basis status is mapped
  to that of the slack. Similarly, the range constraint should be reported infeasible
  if either the slack's bounds or the equality are:

  .. code-block:: c++

    IIS GurobiBackend::GetIIS() {
      pre::ModelValuesInt mv = GetValuePresolver().
        PostsolveIIS( pre::ModelValuesInt{ VarsIIS(), ConsIIS() } );
      return { mv.GetVarValues()(), mv.GetConValues()() };
    }


Implementation API
^^^^^^^^^^^^^^^^^^

To implement value pre- / postsolving, the following API is used:

- `mp::pre::ValuePresolver` implements the interface of
  `mp::pre::BasicValuePresolver`. It calls the individual pre- and postsolve
  routines.

- `mp::pre::BasicLink` is the interface to various implementations of links
  between nodes, such as
  `mp::pre::CopyLink`. Templates `mp::pre::BasicIndivEntryLink` and
  `mp::pre::BasicStaticIndivEntryLink` are base classes for links such as
  `mp::pre::RangeCon2Slack`.


C++ ASL adapter
---------------

An efficient type-safe `C++ adapter for the traditional ASL library
<https://github.com/ampl/mp/tree/master/src/asl>`_ for
connecting solvers to AMPL and other systems. ASL has many additional functions,
such as writing NL files and automatic differentiation.




More details
------------

This section overviews some more details of the API.

For a complete API reference, see the :ref:`index <genindex>`.



Problem representation
~~~~~~~~~~~~~~~~~~~~~~

A standard representation of a model, convenient for intermediate storage.

* `mp::ProblemInfo`, `mp::var::Type`, `mp::obj::Type`, `mp::func::Type`, `mp::ComplInfo`


Expression forest walkers
~~~~~~~~~~~~~~~~~~~~~~~~~

Typesafe expression walkers for models stored in memory.

* `mp::expr::Kind`, `mp::expr::str`, `mp::expr::nl_opcode`

* `mp::BasicExprVisitor`, `mp::ExprVisitor`, `mp::ExprConverter`

* `mp::ProblemFlattener`



Solution status
~~~~~~~~~~~~~~~

* `mp::sol::Status`


Suffixes
~~~~~~~~

* `mp::suf::Kind`, `mp::SuffixDef`

* Standard suffix value enums: `mp::IISStatus`, `mp::BasicStatus`


