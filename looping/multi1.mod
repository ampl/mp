
# ----------------------------------------
# MULTI-COMMODITY FLOW USING
# DANTZIG-WOLFE DECOMPOSITION
# ----------------------------------------

### SUBPROBLEM ###

set ORIG;   # origins
set DEST;   # destinations
set PROD;   # products

param supply {ORIG,PROD} >= 0;  # amounts available at origins
param demand {DEST,PROD} >= 0;  # amounts required at destinations

   check {p in PROD}:
      sum {i in ORIG} supply[i,p] = sum {j in DEST} demand[j,p];

param price_convex;                  # dual price on convexity constr
param price {ORIG,DEST} <= 0.000001; # dual price on shipment limit
param cost {ORIG,DEST,PROD} >= 0;    # shipment costs per unit

var Trans {ORIG,DEST,PROD} >= 0;   # units to be shipped

minimize Artif_Reduced_Cost:
   sum {i in ORIG, j in DEST, p in PROD}
      (- price[i,j]) * Trans[i,j,p] - price_convex;

minimize Reduced_Cost:
   sum {i in ORIG, j in DEST, p in PROD}
      (cost[i,j,p] - price[i,j]) * Trans[i,j,p] - price_convex;

subject to Supply {i in ORIG, p in PROD}:
   sum {j in DEST} Trans[i,j,p] = supply[i,p];

subject to Demand {j in DEST, p in PROD}:
   sum {i in ORIG} Trans[i,j,p] = demand[j,p];

### MASTER PROBLEM ###

param limit {ORIG,DEST} >= 0;   # max shipped on each link

param nPROP integer >= 0;
param prop_ship {ORIG,DEST,1..nPROP} >= -0.000001;
param prop_cost {1..nPROP} >= 0;

            # For each proposal from the subproblem:
            # amount it ships over each link, and its cost

var Weight {1..nPROP} >= 0;
var Excess >= 0;

minimize Artificial: Excess;

minimize Total_Cost:
   sum {k in 1..nPROP} prop_cost[k] * Weight[k];

subject to Multi {i in ORIG, j in DEST}:
   sum {k in 1..nPROP} prop_ship[i,j,k] * Weight[k] - Excess <= limit[i,j];

subject to Convex: sum {k in 1..nPROP} Weight[k] = 1;

### PHASE III PROBLEM ###

param opt_ship {ORIG,DEST} >= -0.000001;

minimize Opt_Cost:
   sum {i in ORIG, j in DEST, p in PROD} cost[i,j,p] * Trans[i,j,p];

subject to Opt_Multi {i in ORIG, j in DEST}:
   sum {p in PROD} Trans[i,j,p] = opt_ship[i,j];
