Backend
=======

Using the `mp::Backend` and the derived classes is now the recommended approach to building
a new solver interface.
It provides a convenient API for common solver options and suffixes.

Backend, MIPBackend
-------------------

* `mp::Backend`, `mp::MIPBackend` standardize some common AMPL app behaviour, such as
  solver messages and status reporting, simplex basis statuses, and suffix I/O


Solver, SolverImpl
------------------

* `mp::Solver` and `mp::SolverImpl` enable very basic standard behaviour (e.g., multiobj, solution output)


Internal API
------------
