var x{1..4} integer >= -2 <= 1;

minimize o: -3*x[1] - 5*x[2] + 2*x[3] + x[4];

s.t. c1: abs(x[1] - x[3]) >= 0;
