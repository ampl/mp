var x in 1..9;
var y >= 1 <= 200;

minimize Obj: x+y;

s.t. ConCtxPos: y>3 ==> (x==2 || x==6
       || x==5
);

## 9 is on boundary.
## See also test_ineq_unify_01.mod
s.t. ConCtxNeg: (x==9 || x==4 || x==7) ==> y >= 4;
