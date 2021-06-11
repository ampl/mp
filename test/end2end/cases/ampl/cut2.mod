
  problem Cutting_Opt;
# ----------------------------------------

param nPAT integer >= 0, default 0;
param roll_width;

set PATTERNS = 1..nPAT;
set WIDTHS;

param orders {WIDTHS} > 0;
param nbr {WIDTHS,PATTERNS} integer >= 0;

   check {j in PATTERNS}: sum {i in WIDTHS} i * nbr[i,j] <= roll_width;

var Cut {PATTERNS} integer >= 0;

minimize Number: sum {j in PATTERNS} Cut[j];

subject to Fill {i in WIDTHS}:
   sum {j in PATTERNS} nbr[i,j] * Cut[j] >= orders[i];


  problem Pattern_Gen;
# ----------------------------------------

param price {WIDTHS} default 0;

var Use {WIDTHS} integer >= 0;

minimize Reduced_Cost:  
   1 - sum {i in WIDTHS} price[i] * Use[i];

subject to Width_Limit:  
   sum {i in WIDTHS} i * Use[i] <= roll_width;
