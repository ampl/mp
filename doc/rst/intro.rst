.. _library-intro:

MP library introduction
=======================

Download and Installation
-------------------------

Downloading solver binaries
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Binaries for the open-source AMPL solvers for major platforms
can be downloaded from AMPL's `Open Source Solvers`__ page.
To use a solver with AMPL, extract the binaries from a downloaded
archive into the AMPL installation directory.

__ https://ampl.com/products/solvers/open-source/

Binaries for AMPL commercial solvers can be downloaded from
`AMPL Portal`__.

__ https://portal.ampl.com/


Building from source
~~~~~~~~~~~~~~~~~~~~

An included `CMake build script`__ can be used to build the MP library,
solver interfaces and function libraries on a wide range of platforms.
You can download CMake for free from https://www.cmake.org/download/.

__ https://github.com/ampl/mp/tree/develop/CMakeLists.txt

CMake workflow
``````````````

CMake works by generating native makefiles or project files that can
be used in the compiler environment of your choice. The typical
workflow starts as follows:

.. code-block:: bash

  git clone https://github.com/ampl/mp.git
  cd mp
  git checkout develop               # For newest code (but possibly work-in-progress)
  git submodule init
  git submodule update
  mkdir build                        # Create a directory to hold the build output
  cd build
  cmake .. -DBUILD=ilogcp,gecode         # Generate native build scripts
      -DCMAKE_OSX_ARCHITECTURES="x86_64" # Select x86_64 architecture on Apple M1,
                                         # required for IBM CPLEX-based drivers
                                         # as of CPLEX Studio 22.1

Note: If the ``arith.h`` file used by default does not match the target architecture,
or if the compiler is not sufficiently compatible with gcc or Microsoft C/C++,
run ``cmake .. -DBUILD=all -DGENERATE_ARITH=true`` to generate an
architecture-specific ``arith.h`` file with ``arithchk``.

If you are on a \*nix system, you should now see a Makefile in the
current directory. Now you can build MP by running ``make``.

Once MP has been built you can invoke ``make test`` to run unit tests.
See also end-to-end testing in section :ref:`howto-test`.

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
a module. See the documentation in :ref:`solver-drivers` and
:file:`solvers/CMakeLists.txt` for module definitions.

By default all modules are disabled and only the main MP libraries are built.
To enable modules, pass their names as a comma-separated list in the ``BUILD``
variable when running CMake::

  cmake .. -DBUILD=ilogcp,gecode

Use ``-DBUILD=all`` to build all modules.

If a module is enabled, its dependencies are automatically downloaded
and built when necessary. For example, enabling the ``gecode`` module
will download the source code of Gecode__ constraint programming solver,
build the solver and its AMPL interface.

__ https://www.gecode.org/

Dependencies of some modules cannot be handled automatically due to
licensing restrictions. If you enable such module, you should have its
dependencies installed on the systems or it will not be built.
For example, if you enable the ``ilogcp`` module, you should have
`IBM ILOG CPLEX Optimization Studio`__ installed.

__ https://www-03.ibm.com/software/products/en/ibmilogcpleoptistud

It is possible to override the automatic detection of dependencies by
adding variables to cmake command line when configuring the project.
The variables to be defined are of the form depname_LIBS
and depname_INCLUDE_DIRS. For example, to build the solvers *copt* or
*xpress* (which do not have automatic dependency detection),
you can use the following::

  cmake .. -DBUILD=copt,xpress
           -DCOPT_LIBS=d:/copt/libs/win64/copt.lib
           -DCOPT_INCLUDE_DIRS=d:/copt/include
           -DXPRESS_INCLUDE_DIRS=../xpress/41.01.01/include
           -DXPRESS_LIBS="../xpress/41.01.01/lib/linux64/libxprs.so.41;../xpress/41.01.01/lib/linux64/libxprl.so.x9.0"

To build the solver *scipmp* statically (on Linux or MacOS), you need to make sure
to add all its dependencies. The following cmake command builds *scipmp* assuming 
that it is installed with the LP solver SoPlex, the arithmetic library GMP and 
the library bliss for symmetry detection (*SCIP* and its dependencies are all 
installed in the standard location)::

  cmake .. -DBUILD=scipmp
           -DSCIP_LIBS="/usr/local/lib/libscip.a;/usr/local/lib/libsoplex.a;/usr/local/lib/libgmp.a;/usr/local/lib/libbliss.a"
           -DSCIP_INCLUDE_DIRS=/usr/local/include

Note: Any other dependency can be used in the same way.

Similarly, for *ortoolsmp* (on MacOS) assuming *ortools* is installed in
the standard location::

  cmake .. -DBUILD=ortoolsmp
           -DORTOOLSMP_INCLUDE_DIRS=/usr/local/include
           -DORTOOLSMP_LIBS=/usr/local/lib/libortools.dylib



Using Eclipse CDT
`````````````````

You can generate `Eclipse CDT <https://www.eclipse.org/cdt/>`_ project files
with CMake::

  cmake -G "Eclipse CDT 4 -  Unix Makefiles"

To get rid of semantic errors reported by Indexer add preprocessor symbols
``_GLIBCPP_USE_NAMESPACES``, ``__GXX_EXPERIMENTAL_CXX0X__`` and ``STAND_ALONE``
in "Project Properties" / "C/C++ Include Files and Symbols" and rebuild
the index.



Building the documentation
``````````````````````````

To build the documentation (automatically, via CMake) you need `Doxygen`__
as well as Python 3.x with Sphinx and Breathe,
see :file:`doc/requirements.txt` (install automatically by
:code:`pip install -r requirements.txt`).
The HTML output is located in :file:`(build folder)/doc/html/index.html`.
To have the alphabetic index automatically generated, install `pandoc`.

__ https://doxygen.nl/

Configure CMake with :code:`-DBUILD_DOC=off` to switch documentation
building off.




FAQ
---



Contributing
------------

Use the ``develop`` branch for new code.

As an example workflow, see :ref:`howto`.


Troubleshooting
---------------

For general questions, email *support /at\\ ampl.com*.

For technical issues submit a ticket at
`https://github.com/ampl/mp/issues <https://github.com/ampl/mp/issues>`_.

Licenses
--------

Copyright (C) 1990 - 2001 Lucent Technologies

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all
copies and that both that the copyright notice and this
permission notice and warranty disclaimer appear in supporting
documentation, and that the name of Lucent or any of its entities
not be used in advertising or publicity pertaining to
distribution of the software without specific, written prior
permission.

LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.


----------------------------------------------------------------------

Copyright (C) 2007 David M. Gay

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that the copyright notice and this permission notice and warranty
disclaimer appear in supporting documentation.

The author disclaims all warranties with regard to this software,
including all implied warranties of merchantability and fitness.
In no event shall the author be liable for any special, indirect or
consequential damages or any damages whatsoever resulting from loss of
use, data or profits, whether in an action of contract, negligence or
other tortious action, arising out of or in connection with the use or
performance of this software.

----------------------------------------------------------------------


Copyright (C) 2022 AMPL Optimization Inc.

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that the copyright notice and this permission notice and warranty
disclaimer appear in supporting documentation.

The author and AMPL Optimization Inc disclaim all warranties with
regard to this software, including all implied warranties of
merchantability and fitness.  In no event shall the author be liable
for any special, indirect or consequential damages or any damages
whatsoever resulting from loss of use	, data or profits, whether in an
action of contract, negligence or other tortious action, arising out
of or in connection with the use or performance of this software.


