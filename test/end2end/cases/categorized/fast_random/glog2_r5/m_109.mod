var x{1..5} integer >= -2 <= 1;

maximize o: 2*x[1] - x[2] - 2*x[3] + 0*x[4] + 5*x[5];

s.t. c1: ((if x[4] > 1 then 1 else 0) + (if x[5] >= 1 then 1 else 0)) >= 1;
s.t. c2: ((if x[4] > 1 then 1 else 0) + (if x[5] >= 1 then 1 else 0)) >= 1;
