var x{1..5} integer >= 0 <= 3;

minimize o: 2*x[1] + x[2] + 5*x[3] - 4*x[4] + 4*x[5];

s.t. c1: not ((x[1] != 4) and (x[4] <= 0));
s.t. c2: abs(x[1] - x[2]) < 3;
