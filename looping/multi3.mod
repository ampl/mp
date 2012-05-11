
# ----------------------------------------
# MULTI-COMMODITY FLOW USING
# DANTZIG-WOLFE DECOMPOSITION
# (one single-product subproblem model)
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

# ----------------------------------------

param prod symbolic within PROD;

var Trans {ORIG,DEST} >= 0;   # units to be shipped

minimize Artif_Reduced_Cost:
   sum {i in ORIG, j in DEST}
      (- price[i,j]) * Trans[i,j] - price_convex[prod];

minimize Reduced_Cost:
   sum {i in ORIG, j in DEST}
      (cost[i,j,prod] - price[i,j]) * Trans[i,j] - price_convex[prod];

subject to Supply {i in ORIG}:
   sum {j in DEST} Trans[i,j] = supply[i,prod];

subject to Demand {j in DEST}:
   sum {i in ORIG} Trans[i,j] = demand[j,prod];


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
