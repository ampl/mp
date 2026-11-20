var x{1..4} integer >= 0 <= 1;

minimize o: 3*x[1] + 4*x[2] - 2*x[3] + 4*x[4];

s.t. c1: ((if x[3] <= 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
