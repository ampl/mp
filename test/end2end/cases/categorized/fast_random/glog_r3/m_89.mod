var x{1..6} integer >= -2 <= 1;

maximize o: x[1] + 4*x[2] - 4*x[3] + 0*x[4] + 0*x[5] + 5*x[6];

s.t. c1: abs(x[6] - x[2]) != 1;
s.t. c2: abs(x[2] - x[4]) == 1;
