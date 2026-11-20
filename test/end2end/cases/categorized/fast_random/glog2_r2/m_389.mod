var x{1..4} integer >= 0 <= 4;

minimize o: 0*x[1] + 6*x[2] + 5*x[3] + 3*x[4];

s.t. c1: abs(x[1] - x[3]) != 1;
s.t. c2: abs(x[1] - x[3]) != 1;
