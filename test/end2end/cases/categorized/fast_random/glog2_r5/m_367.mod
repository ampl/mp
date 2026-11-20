var x{1..5} integer >= -2 <= 2;

minimize o: -3000*x[1] - 333*x[2] - 3000*x[3] - 333*x[4] + 3000*x[5];

s.t. c1: ((if x[4] > 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
