var x{1..3} integer >= -2 <= 1;

minimize o: x[1] - 2*x[2] - 3*x[3];

s.t. c1: abs(x[2] - x[1]) >= 0;
