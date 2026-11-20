var x{1..6} integer >= 0 <= 1;

maximize o: -x[1] + 2*x[2] - 5*x[3] - 2*x[4] - 2*x[5] + 2*x[6];

s.t. c1: abs(x[4] - x[6]) == 0;
