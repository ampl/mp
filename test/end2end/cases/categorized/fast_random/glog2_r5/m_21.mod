var x{1..5} integer >= 0 <= 2;

maximize o: -4*x[1] - 4*x[2] - 3*x[3] + x[4] + 0*x[5];

s.t. c1: abs(x[3] - x[2]) >= 2;
s.t. c2: max(x[1], x[1]) - min(x[1], x[1]) > 2;
s.t. c3: (x[4] >= 2) or (x[3] > 2) or (x[2] > 0);
