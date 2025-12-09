# A basic test for signpow(), with affine/quadratic expressions as arguments.
# Note: does not work for (affine_expr)^2 because it's outmultiplied

var x >=-2 <=3;
var y >=-1.5 <=0.2;
var z >=-1.5 <=0.0025;

maximize Obj:
   2*(x+0.1) - abs(x+0.1)^2.3*(x+0.1) + 13*y - abs(y-0.002)^2*abs(y-0.002)^.3*(y-0.002)^3 + 360*z*y - abs(y)*y^3*y^2*abs(z-0.001)*abs(z-0.001)^4*(z-0.001)*(z-0.001)^4;
