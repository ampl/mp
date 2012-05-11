
# ----------------------------------------
# MULTI-COMMODITY FLOW USING
# DANTZIG-WOLFE DECOMPOSITION
# (one subproblem for each product)
# ----------------------------------------

### SUBPROBLEM ###

set ORIG;   # origins
set DEST;   # destinations
set PROD;   # products

param supply {ORIG,PROD} >= 0;  # amounts available at origins
param demand {DEST,PROD} >= 0;  # amounts required at destinations

   check {p in PROD}:
      sum {i in ORIG} supply[i,p] = sum {j in DEST} demand[j,p];

param price_convex {PROD};           # dual price on convexity constr
param price {ORIG,DEST} <= 0.000001; # dual price on shipment limit
param cost {ORIG,DEST,PROD} >= 0;    # shipment costs per unit

var Trans {ORIG,DEST,PROD} >= 0;   # units to be shipped

minimize Artif_Reduced_Cost {p in PROD}:
   sum {i in ORIG, j in DEST}
      (- price[i,j]) * Trans[i,j,p] - price_convex[p];

minimize Reduced_Cost {p in PROD}:
   sum {i in ORIG, j in DEST}
      (cost[i,j,p] - price[i,j]) * Trans[i,j,p] - price_convex[p];

subject to Supply {i in ORIG, p in PROD}:
   sum {j in DEST} Trans[i,j,p] = supply[i,p];

subject to Demand {j in DEST, p in PROD}:
   sum {i in ORIG} Trans[i,j,p] = demand[j,p];

### MASTER PROBLEM ###

param limit {ORIG,DEST} >= 0;   # max shipped on each link

param nPROP {PROD} integer >= 0;
param prop_ship {ORIG, DEST, p in PROD, 1..nPROP[p]} >= -0.000001;
param prop_cost {p in PROD, 1..nPROP[p]} >= 0;

            # For each proposal from each subproblem:
            # amount it ships over each link, and its cost

var Weight {p in PROD, 1..nPROP[p]} >= 0;
var Excess >= 0;

minimize Artificial: Excess;

minimize Total_Cost:
   sum {p in PROD, k in 1..nPROP[p]} prop_cost[p,k] * Weight[p,k];

subject to Multi {i in ORIG, j in DEST}:
   sum {p in PROD, k in 1..nPROP[p]}
      prop_ship[i,j,p,k] * Weight[p,k] - Excess <= limit[i,j];

subject to Convex {p in PROD}: sum {k in 1..nPROP[p]} Weight[p,k] = 1;
