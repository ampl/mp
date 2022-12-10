var x >=-10, <=20, integer;
var y >=-20, <=30, integer;
var z >=-25, <=49, integer;


minimize Total01:
        10*x + 5*y + 6*z;

subject to Impl01:
        x <= -5 ==> y <= 0 else z >= -15;

subject to OrImpl03:
        x + y >= 2 || (x - y >= -2 ==> y + 2*z <= -15 else x - z == 4);

subject to Sum01:
        x + y + z >= -30;


