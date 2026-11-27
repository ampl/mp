var x{1..6} integer >= 0 <= 4;

maximize o: -4*x[1] - 2*x[2] - 3*x[3] + 3*x[4] + 6*x[5] + 4*x[6];

s.t. c1: abs(x[5] - x[5]) == 1;
s.t. c2: x[3] == 0;
s.t. c3: x[5] == 0;
