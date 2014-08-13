AMPL/MP
=======

.. image:: https://travis-ci.org/ampl/mp.png?branch=master
  :target: https://travis-ci.org/ampl/mp

.. image:: https://ci.appveyor.com/api/projects/status/91jw051om9q8pwt9
  :target: https://ci.appveyor.com/project/vitaut/mp

An open-source library for mathematical programming.

Features
--------

* Reusable high-performance `.nl <https://en.wikipedia.org/wiki/Nl_(format)>`_
  reader

* Efficient type-safe C++ API for connecting solvers to AMPL and other systems:
  `source <https://github.com/ampl/mp/tree/master/src/asl>`_

* Interfaces to solvers supporting
  `AMPL extensions for logic and constraint programming`__:

  __ http://ampl.com/resources/logic-and-constraint-programming-extensions/

  - IBM ILOG CPLEX and CPLEX CP Optimizer
    (`ilogcp <https://github.com/ampl/mp/tree/master/solvers/ilogcp>`_)

  - `Gecode <https://github.com/ampl/mp/tree/master/solvers/gecode>`_

  - `JaCoP <https://github.com/ampl/mp/tree/master/solvers/jacop>`_

* Interfaces to the following solvers:

  - `LocalSolver <https://github.com/ampl/mp/tree/master/solvers/localsolver>`_
  - `Sulum <https://github.com/ampl/mp/tree/master/solvers/sulum>`_

* Interfaces to other solvers via AMPL Solver Library

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

Usage
-----

Binaries for the open-source AMPL solvers and libraries for major platforms
can be downloaded from the AMPL's `Open Source Solvers`__ page.

__ http://ampl.com/products/solvers/open-source/

Building from source
~~~~~~~~~~~~~~~~~~~~

An included `CMake build script`__ can be used to build the MP library,
solver interfaces and function libraries on a wide range of platforms.
You can download CMake for free from
http://www.cmake.org/cmake/resources/software.html.

__ https://github.com/ampl/mp/blob/master/CMakeLists.txt

CMake works by generating native makefiles or project files that can
be used in the compiler environment of your choice. The typical
workflow starts with::

  mkdir build         # Create a directory to hold the build output.
  cd build
  cmake <path/to/mp>  # Generate native build scripts.

where ``<path/to/mp>`` is a path to the ``mp`` repository.

If you are on a *nix system, you should now see a Makefile in the
current directory. Now you can build MP by running ``make``.

Once MP has been built you can invoke ``make test`` to run the tests.

If you use Windows and have Vistual Studio installed, an ``MP.sln`` file
and several ``.vcproj`` files will be created. You can then build them
using Visual Studio or msbuild.

On Mac OS X with Xcode installed, an ``.xcodeproj`` file will be generated.

Optional packages
`````````````````

* To build the amplgsl library you should have the source code for GSL with
  CMake build support in ${MP_DIR}/solvers/amplgsl/gsl. This version of
  GSL is avaialble for download from https://github.com/ampl/gsl.
  You can retrieve MP and GSL sources at the same time using the command

  ::

    git clone --recursive git://github.com/ampl/mp.git

* To build gecode, the AMPL driver for Gecode constraint programming solver,
  you should have Gecode source code in ${AMPL_DIR}/solvers/gecode/lib.
  You can retrieve MP and Gecode sources at the same time using the command

  ::

    git clone --recursive git://github.com/ampl/mp.git

* To build ilogcp, the AMPL driver for IBM ILOG Constraint Programming
  (CP) Optimizer, you should have IBM ILOG CP Optimizer, CPLEX and Concert
  installed. Normally these are installed as parts of IBM ILOG CPLEX
  Optimization Studio. The code has been tested with Optimization Studio
  version 12.4.


Using Eclipse CDT
`````````````````

You can generate Eclipse CDT project files with CMake::

  cmake -G "Eclipse CDT 4 -  Unix Makefiles"

To get rid of semantic errors reported by Indexer add preprocessor symbols
``_GLIBCPP_USE_NAMESPACES``, ``__GXX_EXPERIMENTAL_CXX0X__`` and ``STAND_ALONE``
in "Project Properties" / "C/C++ Include Files and Symbols" and rebuild
the index.

Using Windows SDK
`````````````````

If you want to build MP with the Windows SDK toolchain, use a helper
script run-cmake.bat__ instead of running CMake directly. This script
configures build environment and runs CMake forwarding all command-line
arguments to it, for example::

  support\cmake\run-cmake -G "Visual Studio 10 Win64" .

__ https://github.com/ampl/mp/blob/master/support/cmake/run-cmake.bat

Links
-----
`AMPL home <http://www.ampl.com/>`_ |
`AMPL book <http://ampl.github.io/ampl-book.pdf>`_ |
`Discussion group <https://groups.google.com/group/ampl>`_ |
`SolverStudio for Excel <http://solverstudio.org/languages/ampl/>`_

`AMPL models by Håkan Kjellerstrand <http://www.hakank.org/ampl/>`_
