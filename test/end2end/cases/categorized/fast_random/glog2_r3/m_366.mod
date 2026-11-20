var x{1..3} integer >= -2 <= 2;

minimize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: (if x[2] < -2 then x[2] else x[2]) == -1;
s.t. c2: max(x[3], x[3]) - min(x[2], x[3]) <= -1;
s.t. c3: ((if x[3] != 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
