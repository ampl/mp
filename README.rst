AMPL/MP
=======

.. image:: https://travis-ci.org/ampl/mp.png?branch=master
  :target: https://travis-ci.org/ampl/mp

.. image:: https://ci.appveyor.com/api/projects/status/91jw051om9q8pwt9
  :target: https://ci.appveyor.com/project/vitaut/mp

MP Library is a set of tools recommended to create new AMPL solver interfaces.
`Full documentation. <https://amplmp.readthedocs.io/en/latest/>`__

Features
--------

Reusable building blocks for new interfaces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* High-performance `.nl file reader <https://amplmp.readthedocs.io/en/latest/rst/nl-reader.html>`_
  which is up to `6x faster
  <http://zverovich.net/slides/2015-01-11-ics/socp-reformulation.html#/14>`_
  than the one provided by ASL. It can be used for most efficient translation of NL format into
  solver API.

* Classes `mp::Backend` and `mp::MIPBackend`
  standardize solver behavior such as common options and suffixes
  and are recommended for new interfaces.

* Classes `mp::FlatConverter` and `mp::MIPFlatConverter` facilitate conversion of
  NL expressions which are not natively accepted by a solver into simpler forms.
  `Logical and CP constraints
  <http://ampl.com/resources/logic-and-constraint-programming-extensions/>`__
  are supported.

* Convenience classes `mp::Problem` and `mp::ColProblem` can be used for
  intermediate storage of the NL model.
  `mp::ExprVisitor` and `mp::ExprFlattener` walk NL forest top-down.


Concrete solver interfaces
~~~~~~~~~~~~~~~~~~~~~~~~~~

* Interfaces to solvers with **expression-based APIs.**
  For solvers with an expression-based API,
  NL forests can be efficiently mapped. For example, AMPL expression
  ``exp()`` maps to IBM ILOG Concert's ``IloExponent``. The library
  has the following C++ interfaces of this kind, all of which support
  `AMPL extensions for logic and constraint programming`__:

  __ http://ampl.com/resources/logic-and-constraint-programming-extensions/

  - `Ilogcp <solvers/ilogcp>`_:
    IBM ILOG CPLEX and CPLEX CP Optimizer

  - `Gecode <solvers/gecode>`_

  - `JaCoP <solvers/jacop>`_

  - `LocalSolver <solvers/localsolver>`_

* Interfaces to solvers with **"flat" APIs** (WIP).
  For solvers with more traditional "flat" APIs, class `mp::MIPFlatConverter`
  translates many non-linear AMPL expressions.
  For example, ``max(a, b)`` is translated into a construct meaning
  ``<new var> = max(a, b)``, which is in turn redefined
  into a MIP construct or passed to the solver (Gurobi: `GRBaddgenconstrMax`).
  Currently there are two experimental implementations:

  - `Gurobi <solvers/gurobidirect>`_

  - `IBM ILOG CPLEX <solvers/cplexdirect>`_

Other utilities
~~~~~~~~~~~~~~~

* An efficient type-safe C++ **adapter for the previous ASL library** for connecting solvers to AMPL and other systems:
  `source <src/asl>`_

* `SMPSWriter <solvers/smpswriter>`_,
  a converter from deterministic equivalent of a two-stage stochastic
  programming (SP) problem written in AMPL to an SP problem in SMPS format.

* **End-to-end solver testing script** for testing of various solver features:
  `source <test/end2end>`_


Examples
--------

Reading an .nl file::

  #include "mp/nl.h"
  #include "mp/problem.h"
  
  mp::Problem p;
  ReadNLFile("diet.nl", p);

Usage
-----

Binaries for the open-source AMPL solvers for major platforms
can be downloaded from the AMPL's `Open Source Solvers`__ page.
To use a solver with AMPL, extract the binaries from a downloaded
archive into the AMPL installation directory.

__ http://ampl.com/products/solvers/open-source/

Building from source
~~~~~~~~~~~~~~~~~~~~

An included `CMake build script`__ can be used to build the MP library,
solver interfaces and function libraries on a wide range of platforms.
You can download CMake for free from http://www.cmake.org/download/.

__ CMakeLists.txt

CMake works by generating native makefiles or project files that can
be used in the compiler environment of your choice. The typical
workflow starts with::

  git clone https://github.com/ampl/mp.git
  cd mp
  git submodule init
  git submodule update
  mkdir build  # Create a directory to hold the build output.
  cd build
  cmake -DBUILD=all ..  # Generate native build scripts.

Note: If the ``arith.h`` file used by default does not match the target architecture,
or if the compiler is not sufficiently compatible with gcc or Microsoft C/C++,
run ``cmake -DBUILD=all -DGENERATE_ARITH=true ..`` to generate an architecture-specific ``arith.h`` file with ``arithchk``.

If you are on a \*nix system, you should now see a Makefile in the
current directory. Now you can build MP by running ``make``.

Once MP has been built you can invoke ``make test`` to run the tests.

If you use Windows and have Vistual Studio installed, an ``MP.sln`` file
and several ``.vcproj`` files will be created. You can then build them
using Visual Studio or msbuild.

On Mac OS X with Xcode installed, an ``.xcodeproj`` file will be generated.

Modules
```````

AMPL/MP allows building only parts of the project you are interested in,
for example you can choose to build only a single solver interface.
This is done with the help of modules which are optional components that
can be built separately. Each solver interface and function library is
a module.

By default all modules are disabled and only the main MP libraries are built.
To enable modules, pass their names as a comma-separated list in the ``BUILD``
variable when running CMake::

  cmake -DBUILD=ilogcp,gecode .

Use ``-DBUILD=all`` to build all modules.

If a module is enabled, its dependencies are automatically downloaded
and built when necessary. For example, enabling the ``gecode`` module
will download the source code of Gecode__ constraint programming solver,
build the solver and its AMPL interface.

__ http://www.gecode.org/

Dependencies of some modules cannot be handled automatically due to
licensing restrictions. If you enable such module, you should have its
dependencies installed on the systems or it will not be built.
For example, if you enable the ``ilogcp`` module, you should have
`IBM ILOG CPLEX Optimization Studio`__ installed.

__ http://www-03.ibm.com/software/products/en/ibmilogcpleoptistud


Using Eclipse CDT
`````````````````

You can generate `Eclipse CDT <http://www.eclipse.org/cdt/>`_ project files
with CMake::

  cmake -G "Eclipse CDT 4 -  Unix Makefiles"

To get rid of semantic errors reported by Indexer add preprocessor symbols
``_GLIBCPP_USE_NAMESPACES``, ``__GXX_EXPERIMENTAL_CXX0X__`` and ``STAND_ALONE``
in "Project Properties" / "C/C++ Include Files and Symbols" and rebuild
the index.


Building the documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~

To build the documentation (automatically, via CMake) you need Python 3.7 with Sphinx and Breathe.
The HTML output is located in (build folder)/doc/index.html.


Links
-----
`AMPL home <http://www.ampl.com/>`_ |
`AMPL book <http://ampl.github.io/ampl-book.pdf>`_ |
`Discussion group <https://groups.google.com/group/ampl>`_ |
`SolverStudio for Excel <http://solverstudio.org/languages/ampl/>`_

`AMPL models by Håkan Kjellerstrand <http://www.hakank.org/ampl/>`_
