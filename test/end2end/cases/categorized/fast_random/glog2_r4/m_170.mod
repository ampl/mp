var x{1..6} integer >= -2 <= 2;

minimize o: 2*x[1] - 2*x[2] - 4*x[3] + 6*x[4] + 6*x[5] - x[6];

s.t. c1: ((if x[3] >= 2 then 1 else 0) + (if x[5] >= 2 then 1 else 0)) >= 1;
s.t. c2: not ((x[3] >= -1) and (x[5] <= -2));
