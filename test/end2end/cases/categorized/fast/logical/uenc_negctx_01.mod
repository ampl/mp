var x in 1..9;
var y >= 1 <= 200;

minimize Obj: x+y;

s.t. ConCtxPos: y>3 ==> (x==2 || x==6
       || x==5
);

s.t. ConCtxNeg: x==1 ==> y >= 4;
