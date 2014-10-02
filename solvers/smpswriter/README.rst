smpswriter
==========

smpswriter converts a deterministic equivalent of a two-stage
stochastic programming (SP) problem written in AMPL to an SP problem
in `SMPS format <http://myweb.dal.ca/gassmann/smps2.htm>`__.
It is written as an AMPL solver but can be used as a stand-alone program
to convert `.nl files <http://en.wikipedia.org/wiki/Nl_(format)>`__
(requires .col and .row files as well) into SMPS.

Features
--------

* Supports continuous, binary and integer variables.

* Automatically deduces probabilities from the objective function.

* Supports randomness in

  - constraint right-hand sides
  - constraint matrix elements
  - variable bounds
  
* The same problem can be used for solving a deterministic equivalent in
  AMPL and for generating SMPS using smpswriter.

Problem requirements
--------------------

To be convertible to SMPS with smpswriter the AMPL problem should satisfy
the following requirements:

1. Second-stage variables should be marked with `suffix stage 2`:

   .. code-block:: python

      var sell{Crops, Scenarios} >= 0, suffix stage 2;

2. Second-stage variables and constraints should be indexed over a scenario
   set which should be the last in indexing:

   .. code-block:: python

      set Scenarios;
      var sell{Crops, Scenarios} >= 0, suffix stage 2;
      s.t. quota{s in Scenarios}: sell['beets', s] <= BeetsQuota;

3. Objective should contain at least one second-stage variable, no random
   (scenario-dependent) parameters and expression should have the form of
   expectation:

  .. code-block:: python

    param P{Scenarios}; # probabilities

    maximize profit: sum{s in Scenarios} P[s] * (
      ExcessSellingPrice * sellExcess[s] +
      sum{c in Crops} (SellingPrice[c] * sell[c, s] -
                       PurchasePrice[c] * buy[c, s]) -
      sum{c in Crops} PlantingCost[c] * area[c]);

Usage example
-------------

The file `farmer.ampl <https://raw.github.com/vitaut/ampl/master/solvers/smpswriter/farmer.ampl>`__
contains the deterministic equivalent of the farmer's problem from the
book Introduction to Stochastic Programming by Birge and Louveaux.

1. Read the deterministic equivalent of the farmer's problem written in AMPL
   and write .nl, .col and .row files:

   .. code-block:: python

      ampl: include farmer.ampl;
      ampl: suffix stage IN;    # make the write command output stage information
      ampl: option auxfiles rc; # make the write command output .row and .col files
      ampl: write gfarmer;      # write farmer.nl, farmer.row and farmer.col files

2. Convert the deterministic equivalent in the .nl format into the SP problem
   in SMPS format:

   .. code-block:: shell

      $ smpswriter farmer

   You can do the same from AMPL itself:
   
   .. code-block:: python

      ampl: shell 'smpswriter farmer';

3. Use the generated SMPS files as you like. For this example let's solve the
   SP version of the farmer's problem using `the FortSP solver
   <http://www.optirisk-systems.com/products_fortsp.asp>`__:
   
   .. code-block:: shell

      $ fortsp --smps-obj-sense=maximize farmer
      Stage 1 has 1 row(s), 3 column(s), and 3 nonzero(s).
      Stage 2 has 4 row(s), 7 column(s), and 12 nonzero(s).
      Problem has 2 stage(s) and 3 scenario(s).
      Itn      Objective          Bound        Rel.Gap
        1         107240         136400       0.271914
        2         107240         115429      0.0763595
        3         107240         111053      0.0355573
        4         107240         110011      0.0258422
        5         107240         108861      0.0151167
        6         108328         108802     0.00437892
        7         108328         108390    0.000574388
        8         108390         108390    5.37021e-16
      Number of iterations = 8.
      Master time = 0.001168 s.
      Recourse time = 0.002706 s.
      Optimal solution found, objective = 108390.
      Solution time = 0.011307 s.

   and compare the optimal value to the one found by solving the deterministic
   equivalent:
   
   .. code-block:: text

      ampl: solve;
      MINOS 5.51: optimal solution found.
      11 iterations, objective 108389.8916

Limitations
-----------

* Ranges are not supported.

* Random objective coefficients are not supported. A simple workaround is
  to introduce an auxiliary second-stage variable equal to the objective
  expression (without the expectation) and use this variable in the
  objective function.
