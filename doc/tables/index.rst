Database and spreadsheet connection guide
=========================================

.. epigraph::

   The Guide is definitive. Reality is frequently inaccurate.

   -- Douglas Adams

AMPL's `relational database facility <http://www.ampl.com/NEW/tables.html>`_
provides a convenient mechanism for accessing data in spreadsheets and
databases. A flexible ``table`` declaration defines correspondences between
sets, parameters, variables, and expressions in AMPL models and external
tables managed by database systems or spreadsheet applications.
The ``read table`` and ``write table`` commands use these correspondences to
import data values into AMPL and to export data and solution values from
AMPL.

AMPL can connect to external data sources through a standard `Open Database
Connectivity (ODBC) <http://en.wikipedia.org/wiki/ODBC>`_ interface.
This requires an ODBC driver for the type of data source you are using.
For example if you want to connect to an Oracle database, you should have
an `Oracle ODBC driver
<http://www.oracle.com/technetwork/database/windows/index-098976.html>`_
installed.

There are ODBC drivers for most database management systems (DBMS) and some
non-DBMS data sources such as Microsoft Excel and even CSV files.

This document will guide you through the set up necessary to connect AMPL
to one of the following database systems:

.. toctree::
   :maxdepth: 2

   mysql
