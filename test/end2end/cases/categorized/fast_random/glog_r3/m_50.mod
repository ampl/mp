var x{1..6} integer >= 0 <= 2;

maximize o: -2*x[1] + 2*x[2] - 4*x[3] - 2*x[4] + 2*x[5] + x[6];

s.t. c1: (x[1] <= 0) ==> (x[2] == 1);
s.t. c2: not ((x[6] > 0) and (x[6] > 2));
s.t. c3: max(x[5], x[3]) - min(x[1], x[2]) <= 1;
