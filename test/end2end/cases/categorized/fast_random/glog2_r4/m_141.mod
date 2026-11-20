var x{1..6} integer >= 0 <= 4;

minimize o: 0*x[1] + x[2] - 4*x[3] + 3*x[4] - 3*x[5] + 5*x[6];

s.t. c1: ((if x[3] != -1 then 1 else 0) + (if x[6] >= -1 then 1 else 0)) >= 1;
s.t. c2: (if x[2] == -1 then x[5] else x[2]) > -1;
