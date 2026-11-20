var x{1..3} integer >= 0 <= 3;

minimize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: abs(x[2] - x[3]) != 1;
s.t. c2: not ((x[3] > 1) and (x[3] != 2));
