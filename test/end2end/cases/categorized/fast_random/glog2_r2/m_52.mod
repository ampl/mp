var x{1..5} integer >= 0 <= 2;

minimize o: -3*x[1] - 5*x[2] + 2*x[3] + x[4] - 4*x[5];

s.t. c1: ((if x[1] > 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
