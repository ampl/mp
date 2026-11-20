var x{1..6} integer >= 0 <= 3;

minimize o: -x[1] - 2*x[2] - x[3] + 4*x[4] - 2*x[5] - 2*x[6];

s.t. c1: ((if x[5] <= 3 then 1 else 0) + (if x[1] >= 3 then 1 else 0)) >= 1;
s.t. c2: (if x[2] >= 2 then x[1] else x[3]) == 1;
