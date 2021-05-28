
# -------------------------------------------------------------
# MAXIMUM AND MINIMUM
# maxAndMin1.mod: using <=3 integer expressions
# with "max" and "min" operators
# -------------------------------------------------------------

param n integer := 2;
param ub integer := 10;
param maxValue integer := 8;

set RANGE := 1..n;

var x {RANGE} integer >= 1, <= ub;

maximize TotalCost:
    x[1]-x[2];

subj to MAXIMUM: 
    maxValue = max (x[1]-3, x[2]-2) + 
                    min (x[1]-5, x[2]-6, x[1]-x[2]);
