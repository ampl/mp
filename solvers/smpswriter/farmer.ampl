# The farmer's problem in AMPL
#
# Reference:
#  John R. Birge and Francois Louveaux. Introduction to Stochastic Programming.
#
# AMPL coding by Victor Zverovich.

function expectation;
function random;

suffix stage IN;

set Crops;

set Scen;
param P{Scen}; # probabilities

param TotalArea;               # acre
param PlantingCost{Crops};     # $/acre
param SellingPrice{Crops};     # $/T
param ExcessSellingPrice;      # $/T
param PurchasePrice{Crops};    # $/T
param MinRequirement{Crops};   # T
param BeetsQuota;              # T

# Area in acres devoted to crop c.
var area{c in Crops} >= 0;

# Tons of crop c sold (at favourable price) under scenario s.
var sell{c in Crops} >= 0, suffix stage 2;

# Tons of sugar beets sold in excess of the quota under scenario s.
var sell_excess >= 0, suffix stage 2;

# Tons of crop c bought under scenario s
var buy{c in Crops} >= 0, suffix stage 2;

# The random variable (parameter) representing the yield of crop c.
var RandomYield{c in Crops};

# Realizations of the yield of crop c.
param Yield{c in Crops, s in Scen}; # T/acre

maximize profit:
  expectation(
    ExcessSellingPrice * sell_excess +
    sum{c in Crops} (SellingPrice[c] * sell[c] -
                     PurchasePrice[c] * buy[c])) -
  sum{c in Crops} PlantingCost[c] * area[c];

s.t. totalArea: sum {c in Crops} area[c] <= TotalArea;

s.t. requirement{c in Crops}:
  RandomYield[c] * area[c] - sell[c] + buy[c] >= MinRequirement[c];

s.t. quota: sell['beets'] <= BeetsQuota;

s.t. sellBeets:
  sell['beets'] + sell_excess <= RandomYield['beets'] * area['beets'];

yield: random({c in Crops} (RandomYield[c], {s in Scen} Yield[c, s]));

data;

set Crops := wheat corn beets;
set Scen := below average above;

param TotalArea := 500;

param Yield:
           below average above :=
    wheat    2.0     2.5   3.0
    corn     2.4     3.0   3.6
    beets   16.0    20.0  24.0;

param PlantingCost :=
    wheat 150
    corn  230
    beets 260;

param SellingPrice :=
    wheat 170
    corn  150
    beets  36;

param ExcessSellingPrice := 10;

param PurchasePrice :=
    wheat 238
    corn  210
    beets 100;

param MinRequirement :=
    wheat 200
    corn  240
    beets   0;

param BeetsQuota := 6000;
