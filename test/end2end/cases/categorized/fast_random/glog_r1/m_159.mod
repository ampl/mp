var x{1..4} integer >= -2 <= 2;

minimize o: x[1] + 0*x[2] + x[3] - x[4];

s.t. c1: ((if x[3] < -2 then 1 else 0) + (if x[4] >= -2 then 1 else 0)) >= 1;
