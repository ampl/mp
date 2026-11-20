var x{1..6} integer >= 0 <= 4;

maximize o: 3*x[1] + x[2] + 4*x[3] + 6*x[4] + 5*x[5] - 3*x[6];

s.t. c1: max(x[2], x[1]) - min(x[3], x[2]) == 0;
