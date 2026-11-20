var x{1..4} integer >= -2 <= 2;

minimize o: 6*x[1] + 3*x[2] + 6*x[3] - 3*x[4];

s.t. c1: abs(x[4] - x[1]) <= 1;
