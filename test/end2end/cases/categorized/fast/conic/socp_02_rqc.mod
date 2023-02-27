# socp_02_rqc.mod

var x {1..3} >=0, <=20;

minimize Obj:
   2*x[1] + 3*x[2] - 0.5*x[3];

s.t. RotatedCone:
   25*x[1]*x[2] >= 16*x[3]^2;

s.t. LinCon:
   x[1] + 2*x[2] + 3*x[3] == 6;
