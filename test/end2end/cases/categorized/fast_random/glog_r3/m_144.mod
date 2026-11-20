var x{1..5} integer >= 0 <= 2;

minimize o: -4*x[1] - 5*x[2] + 2*x[3] + 3*x[4] + 5*x[5];

s.t. c1: abs(x[4] - x[2]) == 1;
s.t. c2: abs(x[1] - x[2]) <= 2;
s.t. c3: not ((x[1] >= 2) and (x[5] >= 1));
