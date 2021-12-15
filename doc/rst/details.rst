.. _cppreference:

Details
=======


Common definitions
------------------

Header: :file:`mp/common.h`

Problem information
^^^^^^^^^^^^^^^^^^^

.. doxygenstruct:: mp::ProblemInfo
   :members:

.. doxygenenum:: mp::var::Type

.. doxygenenum:: mp::obj::Type

.. doxygenenum:: mp::func::Type

.. doxygenclass:: mp::ComplInfo
   :members:


Status and suffix values
^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenenum:: mp::sol::Status

.. doxygenenum:: mp::IISStatus

.. doxygenenum:: mp::BasicStatus



Problem builders
----------------

Header: :file:`mp/problem-builder.h`

.. doxygenclass:: mp::ProblemBuilder
   :members:

.. doxygenclass:: mp::ColProblem
   :members:


Header: :file:`mp/problem.h`

.. doxygentypedef:: mp::Problem

.. doxygenclass:: mp::BasicProblem
   :members:


AMPL Expressions
----------------

Header: :file:`mp/basic-expr-visitor.h`

.. doxygenenum:: mp::expr::Kind

.. doxygenfunction:: mp::expr::str

.. doxygenfunction:: mp::expr::nl_opcode

.. doxygenclass:: mp::BasicExprVisitor
   :members:


Header: :file:`mp/expr-visitor.h`

.. doxygenclass:: mp::ExprVisitor
   :members:

.. doxygenclass:: mp::ExprConverter
   :members:


Expression forest flattener
---------------------------

Header: :file:`mp/flat/expr-flattener.h`

.. doxygenclass:: mp::ExprFlattener
   :members:


Solver logic
------------

Header: :file:`mp/solver.h`

.. doxygenclass:: mp::SolverImpl
   :members:

.. doxygenclass:: mp::Solver
   :members:


Suffixes
--------

Header: :file:`mp/suffix.h`


.. doxygenenum:: mp::suf::Kind

.. doxygenclass:: mp::SuffixDef
   :members:


Namespace mp
------------

.. doxygennamespace:: mp
