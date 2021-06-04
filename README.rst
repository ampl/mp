AMPL/MP
=======

.. image:: https://travis-ci.org/ampl/mp.png?branch=master
  :target: https://travis-ci.org/ampl/mp

.. image:: https://ci.appveyor.com/api/projects/status/91jw051om9q8pwt9
  :target: https://ci.appveyor.com/project/vitaut/mp

An open-source library for mathematical programming.


Features
--------

* **Expression-based solver interface.**
  For solvers with an expression API,
  expression trees can be efficiently mapped. For example, AMPL expression
  ``max(a, b)`` directly maps to IBM ILOG Concert's ``IloMax``. The library
  has the following C++ interfaces of this kind, all of which support
  `AMPL extensions for logic and constraint programming`__:

  __ http://ampl.com/resources/logic-and-constraint-programming-extensions/

  - `Ilogcp <solvers/ilogcp>`_:
    IBM ILOG CPLEX and CPLEX CP Optimizer

  - `Gecode <solvers/gecode>`_

  - `JaCoP <solvers/jacop>`_

  - `LocalSolver <solvers/localsolver>`_

* **Conversion-based solver interface (WIP).**
  For solvers with more traditional 'flat' APIs, a customizable conversion
  layer translates expressions into constraints. For example, ``max(a, b)``
  is translated into ``<new var> = max(a, b)``, which is in turn redefined
  into a MIP construct or passed to the solver (Gurobi: ``GRBgenconstrMax``.)
  `Logical and CP constraints
  <http://ampl.com/resources/logic-and-constraint-programming-extensions/>`__
  are supported as well. For the solver API, an easy-to-adapt C++ wrapper class
  is provided. Currently it's two experimental interfaces:
  
  - `Gurobi <solvers/gurobidirect>`_ ('gurobidirect')

  - `IBM ILOG CPLEX <solvers/cplexdirect>`_ ('cplexdirect')

* **End-to-end solver testing script** for testing of various solver features.
  `Documentation. <test/end2end>`_

* An efficient type-safe C++ **adapter for the previous ASL library** for connecting solvers to AMPL and other systems:
  `source <src/asl>`_

* Reusable high-performance `.nl file reader <http://ampl.github.io/nl-reader.html>`_
  which is up to `6x faster
  <http://zverovich.net/slides/2015-01-11-ics/socp-reformulation.html#/14>`_
  than the one provided by ASL. Documentation: http://ampl.github.io

* Database support on Linux and MacOS X.
  See `Database and spreadsheet connection guide`__.

  __  http://ampl.github.io/tables/

* `SMPSWriter <solvers/smpswriter>`_,
  a converter from deterministic equivalent of a two-stage stochastic
  programming (SP) problem written in AMPL to an SP problem in SMPS format.

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

To build the documentation you need the following software installed on your
system:

* `Python <https://www.python.org/>`_ with pip and virtualenv
* `Doxygen <http://www.doxygen.org/>`_

First generate makefiles or project files using CMake as described in
the previous section. Then compile the ``doc`` target/project, for example::

  make doc

This will generate the HTML documenation in ``doc/ampl.github.io``.


Links
-----
`AMPL home <http://www.ampl.com/>`_ |
`AMPL book <http://ampl.github.io/ampl-book.pdf>`_ |
`Discussion group <https://groups.google.com/group/ampl>`_ |
`SolverStudio for Excel <http://solverstudio.org/languages/ampl/>`_

`AMPL models by Håkan Kjellerstrand <http://www.hakank.org/ampl/>`_
