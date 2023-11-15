
MP library
============

MP library is a set of solver drivers and tools recommended to create
new AMPL solver drivers. It provides type-safe and flexible interfaces
suitable for linear and mixed-integer, non-linear, and
Constraint Programming solvers. 

Moreover, MP library provides
a :ref:`lightweight NL writer <write-nl-read-sol>`.

MP replaces the previous
`AMPL Solver Library`__ for solvers not requiring automatic differentiation.
It is implemented in C++, parts of the API are exported in C.

__ https://github.com/ampl/asl


Why MP library?
---------------

MP library aims to provide efficient, simple-to-use and highly
configurable solver drivers and tools. Benefits of MP include:

* Automatic reformulation of some expressions: problems containing such expressions
  can be seamlessly solved using solver drivers developed with MP, although the 
  underlying solver might natively not support them
* Consistency: many solver options are supported via :ref:`features-guide` at framework level, 
  requiring only the implementation of an API in the solver driver. This ensures
  consistency in options naming and semantics across solvers
* Speed of development: modern design patterns and a declarative approach to features
  implementation greatly reduce development effort: for a typical MIP solver, a few days
  are enough to code a driver


On this documentation
----------------------

The documentation is divided in two sections:

* :ref:`userDocumentation` containing the information useful to a modeler/solver driver user.
  
  * :ref:`modeling-guide` presents the advanced reformulation support for many expression types
  * :ref:`features-guide` describes the common features available via solver drivers developed with the framework
  * :ref:`solver-drivers` contains a categorized list of the currently available drivers
* :ref:`developerDocumentation` targeted to developers who want to develop a solver driver 
  
  * :ref:`library-intro` shows instructions for downloading, installing and building the library
  * :ref:`howto` illustrates how to get familiar with the framework and how to quickly obtain the
    boilerplate code necessary to develop a new solver driver
  * :ref:`components` shows the main components of the framework,
    including NL reader and NL writer APIs
  * :ref:`howto-test` illustrates how to execute tests during development
  * :ref:`cppreference` contains the reference to all the classes in the framework 


.. _userDocumentation:

User documentation 
------------------

.. toctree::
   :maxdepth: 2
   
   Modeling guide <model-guide>
   Features guide <features-guide>
   Solver drivers <drivers>

.. _developerDocumentation:

Developer documentation
-----------------------

.. toctree::
   :maxdepth: 3
   
   Developer docs <developers>


Indices and tables
------------------

* :ref:`genindex`
* :ref:`search`
