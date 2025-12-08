# A basic test for signpow()

var x >=-2 <=3;
var y >=-1.5 <=0.2;
var z >=-1.5 <=0.0025;

maximize Obj:
   2*x - abs(x)*x + 13*y - abs(y)^2*abs(y)^.3*y^3 + 360*z*y - abs(y)*y^3*y^2*abs(z)*abs(z)^4*z*z^4;
