## Check that we eliminate top-level AND,
## both for algebraic and logical terms

set I default 1..3;

var x{I} >=0 <=10 integer;

s.t. C1: (x[1]>=3 && x[2]<=5 && x[3]<=7);

s.t. C2: (x[1]<=4 || x[2]>=2) && alldiff(x[2]+7, x[3]+3);

s.t. C3: x[1]+x[2]+x[3] <= 7;

minimize O: 2*x[1]+x[2]-x[3];
