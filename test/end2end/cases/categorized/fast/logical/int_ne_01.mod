
# -------------------------------------------------------------
# x!=y for integer variables
# -------------------------------------------------------------

param ub integer := 10;

var b logical;
var x integer >= -41, <= ub;
var y integer >= -41, <= ub;

maximize TotalSum:
    x+y;

subj to NE:
    x!=y;
