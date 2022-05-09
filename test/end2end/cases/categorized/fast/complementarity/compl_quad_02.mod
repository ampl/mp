### Model compl_quad_01.mod
### Quadratic complementarity, no objective
### Unique solution

var x >= -10, <= 20;
var y >= -15, <= 30;


s.t. LinCon1: y >= x;

s.t. ComplCon1:
    y >= (x-2)^2 complements x >= y^2;
