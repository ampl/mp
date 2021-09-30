AMPL/MP
=======

.. image:: https://travis-ci.org/ampl/mp.png?branch=master
  :target: https://travis-ci.org/ampl/mp

.. image:: https://ci.appveyor.com/api/projects/status/91jw051om9q8pwt9
  :target: https://ci.appveyor.com/project/vitaut/mp

An open-source library for mathematical programming.

Documentation: http://ampl.github.io

Features
--------

* Reusable high-performance `.nl file reader <http://ampl.github.io/nl-reader.html>`_
  which is up to `6x faster
  <http://zverovich.net/slides/2015-01-11-ics/socp-reformulation.html#/14>`_
  than the one provided by ASL

* Efficient type-safe C++ API for connecting solvers to AMPL and other systems:
  `source <https://github.com/ampl/mp/tree/master/src/asl>`_

* Interfaces to solvers supporting
  `AMPL extensions for logic and constraint programming`__:

  __ http://ampl.com/resources/logic-and-constraint-programming-extensions/

  - `Ilogcp <https://github.com/ampl/mp/tree/master/solvers/ilogcp>`_:
    IBM ILOG CPLEX and CPLEX CP Optimizer

  - `Gecode <https://github.com/ampl/mp/tree/master/solvers/gecode>`_

  - `JaCoP <https://github.com/ampl/mp/tree/master/solvers/jacop>`_

  - `LocalSolver <https://github.com/ampl/mp/tree/master/solvers/localsolver>`_

* Interfaces to COIN-OR Solvers (in the `ampl/coin <https://github.com/ampl/coin>`_ repository):

  - Clp
  - Cbc
  - Ipopt
  - Bonmin
  - Couenne

* Interfaces to other solvers via AMPL Solver Library:
  `list of supported solvers <http://ampl.com/products/solvers/all-solvers-for-ampl/>`_

* Cross-platform build support with `CMake <http://www.cmake.org/>`_ and
  continuous integration systems. This includes third-party solvers and
  libraries (COIN-OR solvers with CMake support are available in the
  `ampl/coin <https://github.com/ampl/coin>`_ repository).

* `AMPLGSL <https://github.com/ampl/mp/tree/master/src/gsl>`_, an AMPL function
  library providing access to the GNU Scientific Library (GSL) functions.
  See the `AMPLGSL documentation <http://ampl.github.io/amplgsl>`_.

* Database support on Linux and MacOS X.
  See `Database and spreadsheet connection guide`__.

  __  http://ampl.github.io/tables/

* `SMPSWriter <https://github.com/ampl/mp/tree/master/solvers/smpswriter>`_,
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

__ https://github.com/ampl/mp/blob/master/CMakeLists.txt

CMake works by generating native makefiles or project files that can
be used in the compiler environment of your choice. The typical
workflow starts with::

  mkdir build         # Create a directory to hold the build output.
  cd build
  cmake <path/to/mp>  # Generate native build scripts.

where ``<path/to/mp>`` is a path to the ``mp`` repository.

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

  cmake -DBUILD=gsl,ilogcp .

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


我弄的

