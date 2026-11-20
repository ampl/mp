var x{1..5} integer >= 0 <= 4;

minimize o: -x[1] + x[2] + 5*x[3] + 4*x[4] + 2*x[5];

s.t. c1: abs(x[2] - x[1]) == 1;
s.t. c2: abs(x[1] - x[3]) != 3;
