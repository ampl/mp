var x{1..6} integer >= 0 <= 3;

maximize o: 3*x[1] - 3*x[2] + 2*x[3] + x[4] + 0*x[5] - x[6];

s.t. c1: abs(x[4] - x[3]) <= 2;
s.t. c2: abs(x[4] - x[5]) <= 1;
s.t. c3: ((if x[3] == 2 then 1 else 0) + (if x[6] >= 2 then 1 else 0)) >= 1;
s.t. c4: max(x[2], x[1]) - min(x[3], x[4]) > 2;
