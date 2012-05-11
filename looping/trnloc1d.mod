
# ----------------------------------------
# LOCATION-TRANSPORTATION PROBLEM
# USING BENDERS DECOMPOSITION
# (using dual formulation of subproblem)
# ----------------------------------------

### SUBPROBLEM FOR EXTREME POINT ###

set ORIG;   # shipment origins (warehouses)
set DEST;   # shipment destinations (stores)

param supply {ORIG} > 0;
param demand {DEST} > 0;

param fix_cost {ORIG} > 0;
param var_cost {ORIG,DEST} > 0;

param build {ORIG} binary;  # = 1 iff warehouse built at i

param BIG := 1.0e+9;        # objective > BIG ==> unbounded

var Supply_Price {ORIG} <= 0;
var Demand_Price {DEST};

maximize Dual_Ship_Cost:
   sum {i in ORIG} Supply_Price[i] * supply[i] * build[i] +
   sum {j in DEST} Demand_Price[j] * demand[j];

subj to Unbounded_Check:
   sum {i in ORIG} Supply_Price[i] * supply[i] * build[i] +
   sum {j in DEST} Demand_Price[j] * demand[j] <= BIG;

subj to Dual_Ship {i in ORIG, j in DEST}:
   Supply_Price[i] + Demand_Price[j] <= var_cost[i,j];

### SUBPROBLEM FOR EXTREME RAY ###

maximize Dummy: 0;

subj to Unbounded_Direction:
   sum {i in ORIG} Supply_Price[i] * supply[i] * build[i] +
   sum {j in DEST} Demand_Price[j] * demand[j] = 1;

subj to Feasible_Direction {i in ORIG, j in DEST}:
   Supply_Price[i] + Demand_Price[j] <= 0;

### MASTER PROBLEM ###

param nCUT >= 0 integer;
param cut_type {1..nCUT} symbolic within {"point","ray"};
param supply_price {ORIG,1..nCUT} <= 0.000001;
param demand_price {DEST,1..nCUT};

var Build {ORIG} binary;   # = 1 iff warehouse built at i
var Max_Ship_Cost;

minimize Total_Cost:
   sum {i in ORIG} fix_cost[i] * Build[i] + Max_Ship_Cost;

subj to Cut_Defn {k in 1..nCUT}:
   if cut_type[k] = "point" then Max_Ship_Cost >=
      sum {i in ORIG} supply_price[i,k] * supply[i] * Build[i] +
      sum {j in DEST} demand_price[j,k] * demand[j];
