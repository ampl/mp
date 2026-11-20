var x{1..6} integer >= -2 <= 2;

maximize o: -3*x[1] + 6*x[2] + 3*x[3] - x[4] - 4*x[5] + 5*x[6];

s.t. c1: ((if x[6] != 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
s.t. c2: ((if x[4] >= 1 then 1 else 0) + (if x[5] >= 1 then 1 else 0)) >= 1;
