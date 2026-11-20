var x{1..3} integer >= -2 <= 2;

minimize o: 6*x[1] + x[2] - 2*x[3];

s.t. c1: abs(x[3] - x[1]) > 1;
