var x >= -1 <= 1;
var y >= -1 <= 1;
var z >= -1 <= 1;

minimize Obj: 0.5*x + 0.5*y - 0.9 / (1.0 + exp(-2*x - 1.5*y));

s.t. Con1: x + z == (x^2 + y^2) / (1.0 + exp(-z));
