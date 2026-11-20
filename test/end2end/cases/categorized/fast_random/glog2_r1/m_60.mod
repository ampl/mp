var x{1..5} integer >= 0 <= 3;

minimize o: 0*x[1] - 2*x[2] + 6*x[3] + 5*x[4] + 4*x[5];

s.t. c1: ((if x[3] <= 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
