var x{1..6} integer >= 0 <= 4;

maximize o: -4*x[1] - 4*x[2] - 2*x[3] - 5*x[4] + 0*x[5] + 5*x[6];

s.t. c1: abs(x[3] - x[5]) >= 1;
