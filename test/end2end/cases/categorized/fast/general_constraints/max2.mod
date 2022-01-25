
# -------------------------------------------------------------
# MAXIMUM
# max2.mod: using 2 integer expressions
# with "max" operator
# -------------------------------------------------------------

param n integer := 2;
param ub integer := 10;
param maxValue integer := 5;

set RANGE := 1..n;

var x {RANGE} integer >= 1, <= ub;

maximize TotalCost:
    x[1]-x[2];

subj to MAXIMUM: 
    maxValue = max (x[1]-3, x[2]-2);
