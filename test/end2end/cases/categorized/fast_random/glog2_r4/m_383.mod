var x{1..5} integer >= 0 <= 2;

maximize o: x[1] + 6*x[2] + 0*x[3] - 4*x[4] - 5*x[5];

s.t. c1: (x[3] > 1) or (x[3] != 1) or (x[3] <= 1);
s.t. c2: x[4] == 0;
s.t. c3: x[5] == 2;
