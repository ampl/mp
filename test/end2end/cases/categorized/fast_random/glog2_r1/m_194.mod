var x{1..4} integer >= -2 <= 2;

minimize o: -x[1] - 3*x[2] + 6*x[3] - 3*x[4];

s.t. c1: (if x[4] != -1 then x[2] else x[4]) == -2;
s.t. c2: not ((x[2] == -2) and (x[4] < 1));
