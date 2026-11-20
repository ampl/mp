var x{1..4} integer >= 0 <= 3;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4];

s.t. c1: ((if x[1] > 3 then 1 else 0) + (if x[3] >= 3 then 1 else 0)) >= 1;
