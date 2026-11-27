var x{1..6} integer >= -2 <= 2;

maximize o: x[1] + 2*x[2] + 2*x[3] + 5*x[4] - 4*x[5] + 2*x[6];

s.t. c1: max(x[1], x[2]) - min(x[5], x[5]) > 3;
