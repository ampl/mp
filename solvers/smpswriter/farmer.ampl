# The deterministic equivalent of the farmer's problem from "Introduction
# to Stochastic Programming" by John R. Birge and Francois Louveaux.

set Crops;

set Scenarios;
param P{Scenarios}; # probabilities

param TotalArea;               # acre
param Yield{Crops, Scenarios}; # T/acre
param PlantingCost{Crops};     # $/acre
param SellingPrice{Crops};     # $/T
param ExcessSellingPrice;      # $/T
param PurchasePrice{Crops};    # $/T
param MinRequirement{Crops};   # T
param BeetsQuota;              # T

# Area in acres devoted to crop c
var area{c in Crops} >= 0;

# Tons of crop c sold (at favourable price) under scenario s
var sell{c in Crops, s in Scenarios} >= 0, suffix stage 2;

# Tons of sugar beets sold in excess of the quota under
# scenario s
var sellExcess{s in Scenarios} >= 0, suffix stage 2;

# Tons of crop c bought under scenario s
var buy{c in Crops, s in Scenarios} >= 0, suffix stage 2;

maximize profit: sum{s in Scenarios} P[s] * (
    ExcessSellingPrice * sellExcess[s] +
    sum{c in Crops} (SellingPrice[c] * sell[c, s] -
                     PurchasePrice[c] * buy[c, s]) -
    sum{c in Crops} PlantingCost[c] * area[c]);

s.t. totalArea: sum {c in Crops} area[c] <= TotalArea;

s.t. requirement{c in Crops, s in Scenarios}:
    Yield[c, s] * area[c] - sell[c, s] + buy[c, s]
        >= MinRequirement[c];

s.t. quota{s in Scenarios}: sell['beets', s] <= BeetsQuota;

s.t. sellBeets{s in Scenarios}:
    sell['beets', s] + sellExcess[s]
        <= Yield['beets', s] * area['beets'];

data;

set Crops := wheat corn beets;
set Scenarios := below average above;

param TotalArea := 500;

param P := 
    below   0.333333
    average 0.333333
    above   0.333333;

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
