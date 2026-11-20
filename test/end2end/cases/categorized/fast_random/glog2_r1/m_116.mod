var x{1..4} integer >= 0 <= 2;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4];

s.t. c1: max(x[2], x[4]) - min(x[2], x[4]) == 1;
s.t. c2: abs(x[4] - x[2]) > 1;
s.t. c3: (x[3] > 1) ==> (x[3] > 2);
