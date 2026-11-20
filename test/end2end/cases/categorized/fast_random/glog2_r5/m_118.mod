var x{1..5} integer >= 0 <= 1;

minimize o: x[1] + x[2] + x[3] + x[4] + 6*x[5];

s.t. c1: abs(x[4] - x[1]) < 1;
