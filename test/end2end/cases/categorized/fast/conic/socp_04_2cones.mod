# socp_04_2cones.mod

var x {1..3} >=0, <=20;

minimize Obj:
   2*x[1] + 3*x[2] - 0.5*x[3];

s.t. StdCone:
   sqrt(19.03*x[1]^2 + x[3]^2) <= sqrt(3.72)*x[2];

s.t. RotatedCone:
   -25*x[1]*x[2] <= -16*x[3]^2;

s.t. LinCon:
   sum {i in 1..3} x[i] == 5;
