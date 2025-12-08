# A negative test for signpow()
#Actualy, some are into pow()

var x >=-2 <=3;
var y >=-1.5 <=0.2;
var z >=-1.5 <=0.0025;

maximize Obj:
   2*x - abs(x)*x*x + 13*y - abs(y)^2*abs(y)^.3*y^.3 + 360*z*y - abs(y^.22)*y^3*y^-2*abs(z)*abs(z)^4*z*z^-4
   + 2*x - sqrt(x^2)^2*x + 13*y - sqrt(y*y)^2*sqrt(y*y*y)^.3*y^3 + 360*z*y - abs(y^.33)*y^3*sqrt(y^4.4)*abs(z)*sqrt(z^4)^2*z*z^-4;

s.t. QuickInfeas:
    sqrt(y*y) + sqrt(y*y*y) >= 15;
