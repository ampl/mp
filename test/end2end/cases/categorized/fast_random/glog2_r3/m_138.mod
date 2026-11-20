var x{1..5} integer >= 0 <= 1;

minimize o: -3*x[1] + x[2] - 5*x[3] - 2*x[4] - 4*x[5];

s.t. c1: abs(x[1] - x[3]) <= 0;
s.t. c2: abs(x[4] - x[3]) > 1;
