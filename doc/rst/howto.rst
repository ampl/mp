How to...
=========

The easiest way to getting started developing a new solver driver using ``mp`` is by
looking at the `visitor <https://github.com/ampl/mp/tree/master/solvers/visitor>`_ mock 
driver.

To build it, you can configure the build system of your choiche and specify the cmake variable `BUILD` appropriately::

  mkdir build
  cd build
  cmake .. -DBUILD=visitor
  make

Once built, executing::

  ./visitor modelfilename.nl

will execute the mock driver, which will simply visit the model represented in the nl file.
The visitor source code can be used as a template to create a new driver, as described in the section below.

Create a new driver
-------------------

Either:

1. Copy all the directory visitor into a new directory - and change its name.
2. Rename all occurrences of the word "visitor"

or:

1. Use the file createDriver.py, which does the two items above automatically