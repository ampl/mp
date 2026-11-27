var x{1..3} integer >= -3 <= 1;

minimize o: -2*x[1] + 6*x[2] + 5*x[3];

s.t. c1: not ((x[1] != -2) and (x[2] != -3));
s.t. c2: abs(x[2] - x[3]) == 2;
