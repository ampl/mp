
# -------------------------------------------------------------
# Or, And, Max, Impl, CommonSubExp
# and_or_max_01.mod:
# Two implications reusing a common subexpression
# Nested logical
# -------------------------------------------------------------

param ub integer := 10;

var b logical;
var x integer >= -41, <= ub;
var y integer >= -41, <= ub;

minimize TotalSum:
    b+x -y;

subj to RIMPL: 
    x<=0 ==> b==1;
    
subj to IMPL:
    b==1 ==> y<=0;

subj to logical1:
    (x<=0 or y>=2)  ==>
          (x<=-5 or
              (max(x,y)<=3 and b==0));
