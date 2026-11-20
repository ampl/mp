var x{1..4} integer >= 0 <= 3;

minimize o: -4*x[1] - x[2] + 4*x[3] + 2*x[4];

s.t. c1: max(x[4], x[2]) - min(x[3], x[4]) < 0;
s.t. c2: (if x[2] > 1 then x[2] else x[3]) < 3;
