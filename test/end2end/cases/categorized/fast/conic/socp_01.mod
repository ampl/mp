# socp_01.mod

var x {1..3} >=0, <=20;

minimize Obj:
   2*x[1] + 3*x[2] - 0.5*x[3];

s.t. StdCone:
   -19.03*x[1]^2 - x[3]^2 >= -3.72*x[2]^2;

s.t. LinCon:
   sum {i in 1..3} x[i] == 5;
