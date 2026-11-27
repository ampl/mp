var x{1..6} integer >= 0 <= 3;

maximize o: -3*x[1] + 6*x[2] - x[3] + 5*x[4] + 4*x[5] + 2*x[6];

s.t. c1: max(x[6], x[1]) - min(x[4], x[6]) < 2;
