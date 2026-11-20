var x{1..6} integer >= 0 <= 3;

minimize o: x[1] + x[2] + 0*x[3] + 0*x[4] + 6*x[5] + 6*x[6];

s.t. c1: (if x[5] == 4 then x[4] else x[4]) >= -1;
s.t. c2: not ((x[6] > 4) and (x[1] <= 3));
