Connecting AMPL to MySQL
========================

Installation
------------

GNU/Linux
~~~~~~~~~

Debian-based distributions
``````````````````````````

The following instructions apply to `Debian <http://www.debian.org/>`__
and Debian-based Linux distributions such as `Ubuntu
<http://www.ubuntu.com/>`__ and `Mint <http://linuxmint.com/>`__.

#. Install the MySQL ODBC driver:

   .. code-block:: bash

      $ sudo apt-get install libmyodbc

#. Register the driver:

   .. code-block:: bash

      $ sudo odbcinst -i -d -f /usr/share/libmyodbc/odbcinst.ini

Other distributions
```````````````````

#. Install `unixODBC <http://www.unixodbc.org>`__ following `these instructions
   <http://www.unixodbc.org/download.html>`__.

#. Install and register the MySQL Connector/ODBC following `these instructions
   <http://dev.mysql.com/doc/refman/5.1/en/connector-odbc-installation.html#connector-odbc-installation-binary-unix>`__.

MacOS X
~~~~~~~

#. Install and register the MySQL Connector/ODBC as described `here
   <http://dev.mysql.com/doc/refman/5.1/en/connector-odbc-installation.html#connector-odbc-installation-binary-macosx>`__.

Windows
~~~~~~~

#. Install and register the MySQL Connector/ODBC following `these instructions
   <http://dev.mysql.com/doc/refman/5.1/en/connector-odbc-installation.html#connector-odbc-installation-binary-windows>`__.

Usage
-----

We'll demonstrate usage of MySQL with AMPL on a small example.
For this example we use the diet problem which finds a combination of foods
that satisfies certain nutritional requirements. It is described in
`Chapter 2 of the AMPL book <http://www.ampl.com/BOOK/CHAPTERS/05-tut2.pdf>`__.

We assume that you've installed the MySQL ODBC driver using
the instructions above and have access to the MySQL ``test`` database.

First download the data for the diet problem `diet-mysql.sql
<https://raw.github.com/vitaut/ampl/master/models/tables/diet-mysql.sql>`__
and import it to MySQL:

   .. code-block:: bash

      $ mysql < diet-mysql.sql

Then download the model file `diet.mod
<https://raw.github.com/vitaut/ampl/master/models/tables/diet.mod>`__
and the script file `diet-mysql.run
<https://raw.github.com/vitaut/ampl/master/models/tables/diet-mysql.run>`__.

The script file first reads the model:

   .. code-block:: none

      model diet.mod;

Then it defines a parameter to hold a connection string. Since the connection
parameters are the same for all table declarations in our example this avoids
unnecessary duplication. In this case we specify all the connection parameters
explicitly, instead you can use a DSN file name or ``"DSN=<dsn-name>"``
as a connection string.

   .. code-block:: none

      param ConnectionStr symbolic = "DRIVER=MySQL; DATABASE=test;";

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

TODO: DSN example

SQL statements
--------------

The default `identifier quote character in MySQL
<http://dev.mysql.com/doc/refman/5.0/en/identifiers.html>`__
is the backtick (`````). AMPL's ODBC table handler detects the quote
character automatically and uses it when necessary. This, however, affects
user-supplied SQL statements which are passed to the MySQL ODBC driver as is
and should use the correct quotation. It is possible to enable support for
the ANSI standard quote character (``"``) in MySQL by setting the SQL mode to
`ANSI_QUOTES
<http://dev.mysql.com/doc/refman/5.1/en/server-sql-mode.html#sqlmode_ansi_quotes>`__.

