var x{1..4} integer >= 0 <= 1;

maximize o: 6*x[1] - 5*x[2] + 6*x[3] + 3*x[4];

s.t. c1: ((if x[2] < 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
