var x{1..6} integer >= 0 <= 2;

maximize o: -x[1] + 5*x[2] + x[3] + 6*x[4] - 3*x[5] + 6*x[6];

s.t. c1: max(x[6], x[6]) - min(x[3], x[3]) <= 1;
s.t. c2: x[1] == 0;
