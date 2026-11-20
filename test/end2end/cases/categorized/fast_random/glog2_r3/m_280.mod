var x{1..6} integer >= 0 <= 4;

maximize o: -2*x[1] + 5*x[2] - 4*x[3] + 3*x[4] + 3*x[5] - 2*x[6];

s.t. c1: ((if x[4] > 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
s.t. c2: not ((x[5] < 3) and (x[2] < 3));
