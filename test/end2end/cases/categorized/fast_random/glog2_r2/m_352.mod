var x{1..5} integer >= -2 <= 0;

maximize o: 4*x[1] + 5*x[2] - 5*x[3] + 2*x[4] - 4*x[5];

s.t. c1: (x[1] <= 0) or (x[4] < 0) or (x[2] != -2);
