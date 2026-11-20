var x{1..3} integer >= 0 <= 2;

minimize o: -5*x[1] + 6*x[2] - x[3];

s.t. c1: ((if x[3] < 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
