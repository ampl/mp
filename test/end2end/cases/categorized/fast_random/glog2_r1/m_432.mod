var x{1..3} integer >= 0 <= 4;

minimize o: -x[1] + 2*x[2] + 2*x[3];

s.t. c1: ((if x[2] > 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
