set ORIG;   # origins
set DEST;   # destinations
set PROD;   # products

param supply {ORIG,PROD} >= 0;  # amounts available at origins
param demand {DEST,PROD} >= 0;  # amounts required at destinations

   check {p in PROD}:
      sum {i in ORIG} supply[i,p] = sum {j in DEST} demand[j,p];

param limit {ORIG,DEST} >= 0;   # maximum shipments on routes
param minload > 0;             # minimum nonzero shipment
param maxserve integer > 0;     # maximum destinations served

param vcost {ORIG,DEST,PROD} >= 0; # variable shipment cost on routes
var Trans {i in ORIG,j in DEST,PROD} >= 0 <= limit[i,j];   # units to be shipped

param fcost {ORIG,DEST} >= 0;      # fixed cost for using a route

## Numerically unstable
## (bounds on Ship hard to deduce from later constraints):
## var Ship {i in ORIG, j in DEST} = sum {p in PROD} Trans[i,j,p];

## Define with bounds
var Ship {i in ORIG, j in DEST} >= 0, <= limit[i,j];

s.t. InitShip {i in ORIG, j in DEST}:
   Ship[i,j] = sum {p in PROD} Trans[i,j,p];

minimize Total_Cost:
   sum {i in ORIG, j in DEST, p in PROD} vcost[i,j,p] * Trans[i,j,p]
 + sum {i in ORIG, j in DEST} if Ship[i,j] >= minload then fcost[i,j];

subject to Supply {i in ORIG, p in PROD}:
   sum {j in DEST} Trans[i,j,p] = supply[i,p];

subject to Demand {j in DEST, p in PROD}:
   sum {i in ORIG} Trans[i,j,p] = demand[j,p];

subject to Ship_Range {i in ORIG, j in DEST}:
   Ship[i,j] = 0  or  minload <= Ship[i,j] <= limit[i,j];

subject to Max_Serve {i in ORIG}:
   count {j in DEST} (Ship[i,j] >= minload) <= maxserve;
