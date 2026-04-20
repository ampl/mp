#############################################
## Test inlining of algebraic expressions
## in constraint and objective linear parts
#############################################

var b {1..4} binary;
var x {1..4} in 0..4;

var count2 = count {i in 1..3} (x[i]==i);   # 2 lineq's
var if2 = if b[3] then 2;

minimize Obj:
   x[1] - x[2] + x[3] - x[4]
   + (if b[1] then 5 else 2)                # Operator precedence!
   + count {i in 2..3} (!b[i])
   - count2;

s.t. C1:
   -numberof x[4] in (x[1], x[2], x[3])     # 2 linle, 3 linge
   + if2
   + count {i in 2..3} (!b[i])
   <= 5;

## TODO <=1
s.t. Fix_If2: if2 <= 0;

s.t. Bound_Count2: count2 >= 1;
