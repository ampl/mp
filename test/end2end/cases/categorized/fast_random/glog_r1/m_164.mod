var x{1..6} integer >= 0 <= 2;

minimize o: -2*x[1] - 2*x[2] + 5*x[3] - x[4] + 5*x[5] - 3*x[6];

s.t. c1: ((if x[3] > 2 then 1 else 0) + (if x[4] >= 2 then 1 else 0)) >= 1;
