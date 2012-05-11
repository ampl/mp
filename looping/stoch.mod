
# ----------------------------------------
# STOCHASTIC PROGRAMMING PROBLEM
# USING BENDERS DECOMPOSITION
# ----------------------------------------

### SUBPROBLEM FOR EXTREME POINT ###

set PROD;     # products
param T > 0;  # number of weeks
set SCEN;     # number of scenarios

param rate {PROD} > 0;          # tons per hour produced
param avail {1..T} >= 0;        # hours available in week
param market {PROD,1..T} >= 0;  # limit on tons sold in week

param prodcost {PROD} >= 0;     # cost per ton produced
param invcost {PROD} >= 0;      # carrying cost/ton of inventory

param revenue {PROD,1..T,SCEN} >= 0;  # projected revenue/ton

param prob {SCEN} >= 0, <= 1;
   check: 0.99999 < sum {s in SCEN} prob[s] < 1.00001;

param inv1 {PROD} >= 0;  # inventory at end of first period

var Make {PROD,2..T,SCEN} >= 0;        # tons produced
var Inv {PROD,2..T,SCEN} >= 0;         # tons inventoried
var Sell {p in PROD, t in 2..T, SCEN}  # tons sold
   >= 0, <= market[p,t];

maximize Stage2_Profit:
   sum {s in SCEN} prob[s] *
      sum {p in PROD, t in 2..T} (revenue[p,t,s]*Sell[p,t,s] -
         prodcost[p]*Make[p,t,s] - invcost[p]*Inv[p,t,s]);

subject to Time {t in 2..T, s in SCEN}:
   sum {p in PROD} (1/rate[p]) * Make[p,t,s] <= avail[t];

subject to Balance2 {p in PROD, s in SCEN}:
   Make[p,2,s] + inv1[p] = Sell[p,2,s] + Inv[p,2,s];

subject to Balance {p in PROD, t in 3..T, s in SCEN}:
   Make[p,t,s] + Inv[p,t-1,s] = Sell[p,t,s] + Inv[p,t,s];

### MASTER PROBLEM ###

param inv0 {PROD} >= 0;  # initial inventory

param nCUT >= 0 integer;
param cut_type {1..nCUT} symbolic within {"point","ray"};

param time_price {2..T,SCEN,1..nCUT} >= -0.000001;
param bal2_price {PROD,SCEN,1..nCUT};
param sell_lim_price {PROD,2..T,SCEN,1..nCUT} >= -0.000001;

var Make1 {PROD} >= 0;
var Inv1 {PROD} >= 0;
var Sell1 {p in PROD} >= 0, <= market[p,1];

var Min_Stage2_Profit >= 0;

maximize Expected_Profit:
   sum {s in SCEN} prob[s] *
     sum {p in PROD} (revenue[p,1,s]*Sell1[p] -
        prodcost[p]*Make1[p] - invcost[p]*Inv1[p])
   + Min_Stage2_Profit;

subj to Cut_Defn {k in 1..nCUT}:
   if cut_type[k] = "point" then Min_Stage2_Profit <=
      sum {t in 2..T, s in SCEN} time_price[t,s,k] * avail[t] +
      sum {p in PROD, s in SCEN} bal2_price[p,s,k] * (-Inv1[p]) +
      sum {p in PROD, t in 2..T, s in SCEN}
         sell_lim_price[p,t,s,k] * market[p,t];

subject to Time1:
   sum {p in PROD} (1/rate[p]) * Make1[p] <= avail[1];

subject to Balance1 {p in PROD}:
   Make1[p] + inv0[p] = Sell1[p] + Inv1[p];

