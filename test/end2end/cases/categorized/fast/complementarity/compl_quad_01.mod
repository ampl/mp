### Model compl_quad_01.mod
### Quadratic complementarity with linear objective

var x >= -10, <= 20;
var y >= -15, <= 30;

maximize LinObj: x+y;

s.t. LinCon1: y <= 1/3.0*(x+1);

s.t. LinCon2: x <= -1/3.0*(y-2);

s.t. ComplCon1:
    y >= -(x+1)^2 - 1 complements y >= x^2 - 1;
