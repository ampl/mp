var x{1..6} integer >= 0 <= 4;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5] + 0*x[6];

s.t. c1: abs(x[3] - x[5]) != 0;
s.t. c2: abs(x[6] - x[5]) >= 2;
