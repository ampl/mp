.. index:: Oracle

Connecting AMPL to Oracle
=========================

To use Oracle with AMPL, you need to have the `Oracle ODBC driver
<http://www.oracle.com/technetwork/database/windows/index-098976.html>`__
installed and to have access to a database server, which could be either
local or remote.

Installation
------------

GNU/Linux
~~~~~~~~~

Debian-based distributions
``````````````````````````

The following instructions apply to `Debian <http://www.debian.org/>`__
and Debian-based Linux distributions such as `Ubuntu
<http://www.ubuntu.com/>`__ and `Mint <http://linuxmint.com/>`__.

#. Download RPM packages of Oracle Instant Client for your Linux platform from the
   `Instant Client Downloads page
   <http://www.oracle.com/technetwork/database/features/instant-client/index-097480.html>`__.
   You will need Basic and ODBC packages. In the Usage example we also use ``sqlplus`` from
   the SQL*Plus package.

#. Convert packages from RPM to DEB format:

   .. code-block:: bash

      sudo alien oracle-instantclient*.rpm
   
#. Install the packages:

   .. code-block:: bash

      sudo dpkg -i oracle-instantclient*.deb

#. Install dependencies:

   .. code-block:: bash

      sudo apt-get install libaio1 unixodbc

#. Add the Oracle library directory to the library search paths:

   .. code-block:: bash

      echo /usr/lib/oracle/*/client*/lib | sudo tee -a /etc/ld.so.conf
      sudo ldconfig

#. Register the ODBC driver:

   .. code-block:: bash

      sudo /usr/share/oracle/*/client*/odbc_update_ini.sh / /usr/lib/oracle/*/client*/lib

#. Create a symbolic link for ``libodbcinst.so`` on ``x86-64``:

   .. code-block:: bash

      sudo ln -s /usr/lib/x86_64-linux-gnu/libodbcinst.so /usr/lib/x86_64-linux-gnu/libodbcinst.so.2

   Use the following command with the ``x86`` version instead:

   .. code-block:: bash

      sudo ln -s /usr/lib/libodbcinst.so /usr/lib/libodbcinst.so.2

#. Set the ``ORACLE_HOME`` environment variable:

   .. code-block:: bash

      export ORACLE_HOME=<installation-dir>

   replacing ``<installation-dir>`` with the actual installation directory which can
   printed with the command ``echo /usr/lib/oracle/*/client*``.

   Alternatively you can add the line ``ORACLE_HOME=<installation-dir>`` to
   ``~/.pam_environment`` to set this environment variable permanently for
   the current user. Use ``/etc/environment`` instead of ``~/.pam_environment``
   for system-wide environment variables.
   See also `Persistent environment variables
   <https://help.ubuntu.com/community/EnvironmentVariables#Persistent_environment_variables>`__.

Go to :ref:`oracle-usage`.

..
  Other distributions
  ```````````````````

  #. Install `unixODBC <http://www.unixodbc.org>`__ following `these instructions
    <http://www.unixodbc.org/download.html>`__.

  #. Install the MySQL Connector/ODBC following `these instructions
    <http://dev.mysql.com/doc/refman/5.1/en/connector-odbc-installation.html#connector-odbc-installation-binary-unix>`__.
    Make sure that you use compatible versions of the ODBC driver
    (Connector/ODBC) and the MySQL client library, otherwise the driver
    library will not load and any connection attempt will fail.

  #. Register the ODBC driver. The easiest way to register the driver is
    by using the ``myodbc-installer`` utility included in the distribution,
    for example:

    .. code-block:: bash

	$ sudo myodbc-installer -d -a -n "MySQL" \
	    -t "DRIVER=/usr/local/lib/libmyodbc5a.so"

    ``/usr/local/lib/libmyodbc5a.so`` is the path to the driver library
    that you installed in the previous step. You might need to change it
    if you have a different version of the driver or installed it in a
    different location. See the name of the ``.so`` file in the ``lib``
    directory of the installation package.

    Note that the MySQL ODBC/Connector distribution does not include a
    setup library. If you invoke ``myodbc-installer --help``, you may see an
    outdated example with a ``SETUP`` attribute specifying a setup library.
    Omit this attribute during the driver registration unless you have
    installed a setup library from some other source.

  Go to :ref:`oracle-usage`.

  MacOS X
  ~~~~~~~

  The easiest way to install the MySQL ODBC driver on Mac is by using an
  installer available for download as a DMG archive from the
  `Connector/ODBC download page on the MySQL website
  <http://dev.mysql.com/downloads/connector/odbc/#downloads>`__.

  Alternatively you can install the MySQL Connector/ODBC as described `here
  <http://dev.mysql.com/doc/refman/5.1/en/connector-odbc-installation.html#connector-odbc-installation-binary-macosx>`__,
  skipping the outdated last step (driver registration) and then register
  the driver with the following command:

  .. code-block:: bash

    $ sudo myodbc-installer -d -a -n "MySQL" \
	-t "DRIVER=/usr/local/lib/libmyodbc5w.so"

  ``/usr/local/lib/libmyodbc5w.so`` is the path to the driver library
  that you installed in the previous step. You might need to change it
  if you have a different version of the driver or installed it in a
  different location. See the name of the ``.so`` file in the ``lib``
  directory of the installation package.

  Note that the MySQL ODBC/Connector distribution does not include a
  setup library. If you invoke ``myodbc-installer --help``, you may see an
  outdated example with a ``SETUP`` attribute specifying a setup library.
  Omit this attribute during the driver registration unless you have
  installed a setup library from some other source.

  Go to :ref:`oracle-usage`.

Windows
~~~~~~~

The ODBC driver for Oracle often comes installed by default on modern versions
of Windows. You can check if the driver is installed by running the ODBC Data Source
Administrator, ``odbcad32.exe``, and looking for Oracle in the ``Drivers`` tab.

.. image:: ../img/odbcad32-oracle.png

If the driver is missing, download one from `OracleODBC Drivers Download Page
<http://www.oracle.com/technetwork/database/windows/downloads/index-096177.html>`__
and install it.

.. _oracle-usage:

Usage
-----

We'll demonstrate usage of Oracle with AMPL on a small example.
For this example we use the diet problem, which finds a combination of foods
that satisfies certain nutritional requirements. It is described in
`Chapter 2 of the AMPL book <http://www.ampl.com/BOOK/CHAPTERS/05-tut2.pdf>`__.

We assume that you've already installed the Oracle ODBC driver using
the instructions above and have access to a local Oracle database.

First download the data for the diet problem `diet-oracle.sql
<http://ampl.github.io/models/tables/diet-oracle.sql>`__
and import it into an Oracle database:

.. code-block:: bash

   $ sqlplus <username>/<password> @diet-oracle.sql

where ``<username>`` is the name of a database user and ``<password>`` is the
user's password.
 
Then download the model file `diet.mod
<http://ampl.github.io/models/tables/diet.mod>`__
and the script file `diet-oracle.run
<http://ampl.github.io/models/tables/diet-oracle.run>`__.

The script file first reads the model:

.. code-block:: none

   model diet.mod;

Then it defines a parameter to hold a connection string. Since the connection
parameters are the same for all table declarations in our example, we
avoid unnecessary duplication. In this case we specify all the connection
parameters explicitly. Alternatively, you could use a DSN file name or
``"DSN=<dsn-name>"`` as a connection string.

.. code-block:: none

   param ConnectionStr symbolic = "DRIVER=Oracle; SERVER=localhost;";

If you are using Linux or MacOS X and have chosen a driver name other
than ``Oracle``, you will have to specify this name instead of ``Oracle``
in the ``DRIVER=Oracle`` attribute in the connection string.

A driver name is chosen automatically during installation on Windows,
so if you are using this OS, you will have to find the driver name and
specify it instead of ``Oracle`` in the connection string.
To discover the driver name on Windows, run the ODBC Data Source
Administrator, ``odbcad32.exe``.  Go to the ``Drivers`` tab where all the
installed drivers are listed and look for the one containing ``Oracle``:

.. image:: ../img/odbcad32-oracle.png

A driver name containing a semicolon (``;``) should be surrounded with
``{`` and ``}`` in a connection string, for example:

.. code-block:: none

   param ConnectionStr symbolic =
     "DRIVER={Oracle ODBC Driver; version 6.01}; SERVER=localhost;";

Next there are several table declarations that use the ``ConnectionStr``
parameter defined previously:

.. code-block:: none

   table dietFoods "ODBC" (ConnectionStr) "Foods":
       FOOD <- [FOOD], cost IN, f_min IN, f_max IN,
       Buy OUT, Buy.rc ~ BuyRC OUT, {j in FOOD} Buy[j]/f_max[j] ~ BuyFrac;

   table dietNutrs IN "ODBC" (ConnectionStr) "Nutrients": NUTR <- [NUTR], n_min, n_max;
   table dietAmts IN "ODBC" (ConnectionStr) "Amounts": [NUTR, FOOD], amt;

Finally the script reads the data from the tables

.. code-block:: none

   read table dietFoods;
   read table dietNutrs;
   read table dietAmts;

solves the problem
                  
.. code-block:: none

   solve;

and writes the solution back to the database:

.. code-block:: none

   write table dietFoods;

Note that the same table ``dietFoods`` is used both for input and output.

Running the ``diet-oracle.run`` script with ampl shows that data connection
is working properly and the problem is easily solved:

.. code-block:: bash

   $ ampl diet-oracle.run
   MINOS 5.51: optimal solution found.
   13 iterations, objective 118.0594032

..
  You can use various database tools such as `MySQL workbench
  <https://www.mysql.com/products/workbench/>`__ or `MySQL command-line tool
  <http://dev.mysql.com/doc/refman/5.5/en/mysql.html>`__ to view the data
  exported to the database from the AMPL script:

  .. image:: ../img/mysql-workbench.png

  SQL statements
  --------------

  The default `identifier quote character in MySQL
  <http://dev.mysql.com/doc/refman/5.0/en/identifiers.html>`__
  is the backquote (`````). AMPL's ODBC table handler detects the quote
  character automatically and uses it when necessary. However,
  user-supplied SQL statements are passed to the MySQL ODBC driver as is
  and should use the correct quotation. It is possible to enable support for
  the ANSI standard quote character (``"``) in MySQL by setting the SQL mode to
  `ANSI_QUOTES
  <http://dev.mysql.com/doc/refman/5.1/en/server-sql-mode.html#sqlmode_ansi_quotes>`__.

  Example:

    .. code-block:: none

	table Foods "ODBC" "DRIVER=MySQL; DATABASE=test;"
	  "SQL=SELECT `FOOD`, `cost` FROM `Foods`;": [FOOD], cost;

Troubleshooting
---------------

This section lists common problems with possible solutions.

The first thing to do in case of an error is to get additional information.
Add the option ``"verbose"`` to the table declaration that causes the error,
for example:

.. code-block:: none

   table dietFoods "ODBC" (ConnectionStr) "Foods" "verbose":
     ...

Then rerun your code and you should get a more detailed error message.

Data source name not found
~~~~~~~~~~~~~~~~~~~~~~~~~~

Verbose error:

.. code-block:: none

   SQLDriverConnect returned -1
   sqlstate = "IM002"
   errmsg = "[unixODBC][Driver Manager]Data source name not found, and no default driver specified"
   native_errno = 0

If the data source name (DSN) was not found as in the example above check 
if it is spelled correctly in the connection string. If you are not using a
DSN, check the driver name instead. On a Unix-based system you can get the
list of installed ODBC drivers using  the following commands:

.. code-block:: none

   $ odbcinst -d -q

On Windows use the ODBC Data Source Administrator (see :ref:`oracle-usage`).

If the driver name contains a semicolon (``;``), check that the name is
surrounded with ``{`` and ``}`` in the connection string, for example:

.. code-block:: none

   table Foods "ODBC" "DRIVER={Oracle ODBC Driver; version 6.01}; DATABASE=test;":
     ...

Driver's SQLAllocHandle on SQL_HANDLE_HENV failed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Verbose error:

.. code-block:: none

   SQLDriverConnect returned -1
   sqlstate = "IM004"
   errmsg = "[unixODBC][Driver Manager]Driver's SQLAllocHandle on SQL_HANDLE_HENV failed"
   native_errno = 0

This error may occur if the ``ORACLE_HOME`` environment variable is not set.
