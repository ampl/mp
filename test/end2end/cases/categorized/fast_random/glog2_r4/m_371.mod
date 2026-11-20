var x{1..4} integer >= 0 <= 4;

minimize o: 5*x[1] + 4*x[2] + 4*x[3] - 2*x[4];

s.t. c1: ((if x[4] == 5 then 1 else 0) + (if x[3] >= 5 then 1 else 0)) >= 1;
