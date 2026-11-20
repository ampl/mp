var x{1..4} integer >= -2 <= 0;

maximize o: 2000*x[1] + 51*x[2] - 333*x[3] + 2000*x[4];

s.t. c1: ((if x[1] >= -2 then 1 else 0) + (if x[2] >= -2 then 1 else 0)) >= 1;
s.t. c2: abs(x[2] - x[3]) <= 1;
