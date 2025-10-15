var x >=0 <=10 integer;
var y >=0 <=11 integer;
var z >=-3 <=7;

s.t. LC1:
   x<=0
   || x-3*z>=19-1e-10
   || -2*y+2*z==-28+1e-5
   || !(x+5*y<65)
   || y+z<=2.4    ## not on bound
   || z==3.8;     ## not on bound

s.t. C1: x+y+z >= 15;

minimize OO: x-y-z;
