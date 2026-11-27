var x{1..4} integer >= 0 <= 3;

maximize o: -5*x[1] - 4*x[2] + x[3] + 6*x[4];

s.t. c1: ((if x[4] >= 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
