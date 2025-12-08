# A similar basic test for signpow() from sqrt()

var x >=-2 <=3;
var y >=-1.5 <=0.2;
var z >=-1.5 <=0.0025;

maximize Obj:
   2*x - sqrt(x^2)*x + 13*y - sqrt(y*y)^2*sqrt(y*y)^.3*y^3 + 360*z*y - abs(y)*y^3*sqrt(y^4)*abs(z)*sqrt(z^4)^2*z*z^4;
