Reference
=========

This section overviews most low-level classes and functions in namespace :cpp:any:`mp`.



Problem builders
----------------

* `mp::ProblemBuilder`

* `mp::ColProblem`


Problem representation
----------------------

A standard representation of a model, convenient for intermediate storage.
Can be converted into solver API by a subclassed `mp::ExprVisitor`.

* `mp::ProblemInfo`, `mp::var::Type`, `mp::obj::Type`, `mp::func::Type`, `mp::ComplInfo`

* `mp::Problem`, `mp::BasicProblem`


Expression forest walkers
-------------------------

Typesafe expression walkers for models stored in memory.

* `mp::expr::Kind`, `mp::expr::str`, `mp::expr::nl_opcode`

* `mp::BasicExprVisitor`, `mp::ExprVisitor`, `mp::ExprConverter`

* `mp::ExprFlattener`


Standard AMPL solver logic
--------------------------

* `mp::Solver`, `mp::SolverImpl`

* `mp::Backend`, `mp::MIPBackend`


Solution status
---------------

* `mp::sol::Status`


Suffixes
--------

* `mp::suf::Kind`, `mp::SuffixDef`

* Standard suffix value enums: `mp::IISStatus`, `mp::BasicStatus`


