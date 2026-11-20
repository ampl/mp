var x{1..6} integer >= 0 <= 4;

maximize o: -2*x[1] + 2*x[2] + 6*x[3] + 5*x[4] + 2*x[5] + 4*x[6];

s.t. c1: abs(x[1] - x[3]) != 3;
