var x{1..4} integer >= 0 <= 3;

minimize o: -3*x[1] - 5*x[2] - x[3] + 0*x[4];

s.t. c1: abs(x[4] - x[1]) != 3;
s.t. c2: abs(x[4] - x[1]) != 3;
