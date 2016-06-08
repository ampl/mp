smpswriter
==========

smpswriter converts a two-stage stochastic programming (SP) problem from AMPL
to `SMPS format <http://myweb.dal.ca/gassmann/smps2.htm>`__.
It is written as an AMPL solver but can be used as a stand-alone program
to convert `.nl files <http://en.wikipedia.org/wiki/Nl_(format)>`__ into SMPS.

Features
--------

* Supports continuous, binary and integer variables.

* Supports randomness in

  - constraint right-hand sides
  - constraint matrix elements

Problem requirements
--------------------

To be convertible to SMPS with smpswriter the AMPL problem should satisfy
the following requirements:

1. Second-stage variables should be marked with `suffix stage 2`:

   .. code-block:: python

      var sell{Crops} >= 0, suffix stage 2;

2. Random variables/parameters should be specified using the random function:

   .. code-block:: python

      # The random variable (parameter) representing the yield of crop c.
      var RandomYield{c in Crops};

      # Realizations of the yield of crop c.
      param Yield{c in Crops, s in Scen}; # T/acre

      yield: random({c in Crops} (RandomYield[c], {s in Scen} Yield[c, s]));

3. Objective expression should have the form of expectation:

   .. code-block:: python

      maximize profit:
        expectation(
          ExcessSellingPrice * sell_excess +
          sum{c in Crops} (SellingPrice[c] * sell[c] -
                           PurchasePrice[c] * buy[c])) -
        sum{c in Crops} PlantingCost[c] * area[c];

Usage example
-------------

The file `farmer.ampl <https://raw.github.com/vitaut/ampl/master/solvers/smpswriter/farmer.ampl>`__
contains the farmer's problem from the book Introduction to Stochastic
Programming by Birge and Louveaux.

1. Read the deterministic equivalent of the farmer's problem written in AMPL
   and write the .nl file:

   .. code-block:: python

      ampl: include farmer.ampl;
      ampl: write gfarmer; # write farmer.nl

2. Convert the problem from .nl to SMPS format:

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

Limitations
-----------

The following constructs are not supported:

* Ranges
* Random bounds
* Random objective coefficients

Randomness in the objective can be represented by introducing an auxiliary
second-stage variable equal to the objective expression (without the
expectation) and using this variable in the objective function.
