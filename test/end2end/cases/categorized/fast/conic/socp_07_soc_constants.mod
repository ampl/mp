# socp_07_soc_constants.mod

var x {1..3} >=0, <=20;

minimize Obj:
   2*x[1] + 3*x[2] - 0.5*x[3];

s.t. StdCone1:
   -19.03*x[1]^2 - x[3]^2 >= -3.72*x[2]^2;

s.t. StdCone2:
   -19.03*x[1]^2 - x[3]^2 >= -3.72;

s.t. StdCone3:
   19.03 + x[3]^2 <= 3.72*x[2]^2;

s.t. StdCone4:
   0 >= -3.72*x[2]^2;

s.t. LinCon:
   sum {i in 1..3} x[i] == 5;
