# socp_05_rqcSqrt_abs.mod

var x {i in 1..3}
   >=if i<=2 then 0 else -40,
   <=20;

minimize Obj:
   2*x[1] + 3*x[2] - 0.5*x[3];

s.t. StdCone:
   sqrt(19.03*x[1]^2 + x[3]^2) <= sqrt(3.72)*x[2];

s.t. RotatedCone:
   5*sqrt(2.7*x[1]*x[2]) >= abs(x[3]);

s.t. LinCon:
   sum {i in 1..3} x[i] == 5;
