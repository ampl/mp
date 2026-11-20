var x{1..3} integer >= 0 <= 4;

minimize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: ((if x[1] < 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
