# socp_08_sqrt_constants.mod

var x {i in 1..3}
   >=if i<=2 then 0 else -40,
   <=20;

minimize Obj:
   2*x[1] + 3*x[2] - 0.5*x[3];

s.t. StdCone1:
   sqrt(19.03*x[1]^2 + x[3]^2) <= sqrt(3.72)*x[2];

s.t. StdCone2:
   -sqrt(19.03 + x[3]^2) >= -sqrt(3.72)*x[2];

s.t. RotatedCone1:
   5*sqrt(1.7*x[1]*x[2]) >= sqrt(x[3]^2 + 4.8);

s.t. RotatedCone2:
   1.45*sqrt(0.7*x[2]) >= sqrt(x[3]^2 + 4.8);


s.t. LinCon:
   sum {i in 1..3} x[i] == 5;
