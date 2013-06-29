# Portfolio selection model with SSD constraints.

include ssd.ampl;

param NumScenarios;
param NumAssets;

set Scenarios = 1..NumScenarios;
set Assets = 1..NumAssets;

# Return of asset a in senario s.
param Returns{a in Assets, s in Scenarios};

# Reference return in scenario s.
param Reference{s in Scenarios};

# Fraction of the budget to invest in asset a.
var invest{a in Assets} >= 0 <= 1 := 1 / NumAssets;

subject to ssd_constraint{s in Scenarios}:
  ssd_uniform(sum{a in Assets} Returns[a, s] * invest[a], Reference[s]);

subject to budget: sum{a in Assets} invest[a] = 1;

data ssd-portfolio-data.ampl;

option solver ssdsolver;
option ssdsolver_options 'scaled=1 abs_tolerance=1e-7';
solve;
display {a in Assets: invest[a] > 0} invest[a];