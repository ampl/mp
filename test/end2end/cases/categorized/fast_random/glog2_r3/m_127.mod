var x{1..4} integer >= 0 <= 2;

maximize o: -3000*x[1] - 999*x[2] - 999*x[3] + 17*x[4];

s.t. c1: (x[1] < 2) or (x[2] > 2) or (x[3] <= 2);
