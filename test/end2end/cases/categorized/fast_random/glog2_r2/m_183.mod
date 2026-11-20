var x{1..4} integer >= 0 <= 4;

minimize o: x[1] + 3*x[2] + x[3] + 3*x[4];

s.t. c1: ((if x[4] == 4 then 1 else 0) + (if x[4] >= 4 then 1 else 0)) >= 1;
s.t. c2: x[2] == 2;
