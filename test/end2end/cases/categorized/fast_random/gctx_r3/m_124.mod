var x{1..5} integer >= 0 <= 2;

maximize o: -3*x[1] - 2*x[2] + 4*x[3] + 4*x[4] - 2*x[5];

s.t. c1: (min(x[1], x[3]) > 0) ==> ((x[5] >= 1) or (x[2] + x[5] < 1));
s.t. c2: not (x[5] >= 1);
