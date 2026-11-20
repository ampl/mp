var x{1..6} integer >= 0 <= 2;

minimize o: -2*x[1] - 3*x[2] + x[3] - x[4] - 4*x[5] + x[6];

s.t. c1: abs(x[1] - x[2]) <= 0;
s.t. c2: (if x[3] <= 1 then x[1] else x[1]) == 2;
s.t. c3: max(x[4], x[5]) - min(x[4], x[3]) != 0;
