## Test inequalities of int vars to non-int values

var b{1..4}: binary;

var x >=3, <=6;

minimize Obj:
    if b[1]>0.5 then 0 else 5*x
    + if b[2]<0.5 then 0 else 5*x
    + if b[3]>=0.5 then 5*x
    + if b[4]<=0.5 then 5*x;
