
# ----------------------------------------
# LOCATION-TRANSPORTATION PROBLEM
# USING LAGRANGIAN RELAXATION
# ----------------------------------------

set CITY;

param build_limit integer;

param demand {i in CITY} integer > 0;
param supply {CITY} integer > 0;

param ship_cost {i in CITY, j in CITY} >= 0;

param mult {CITY} >= 0;   # Lagrange multipliers for Demand constr

var Build {CITY} integer >= 0 <= 1;  # = 1 iff warehouse built at i
var Ship {CITY,CITY} >= 0;           # amounts shipped

minimize Lagrangian:
   sum {i in CITY, j in CITY} ship_cost[i,j] * Ship[i,j] +
   sum {j in CITY} mult[j] * (demand[j] - sum {i in CITY} Ship[i,j]);

minimize Shipping_Cost:
   sum {i in CITY, j in CITY} ship_cost[i,j] * Ship[i,j];

subj to Supply {i in CITY}:
   sum {j in CITY} Ship[i,j] <= supply[i] * Build[i];

subj to Demand {j in CITY}:
   sum {i in CITY} Ship[i,j] >= demand[j];

subj to Limit:  sum {i in CITY} Build[i] <= build_limit;
