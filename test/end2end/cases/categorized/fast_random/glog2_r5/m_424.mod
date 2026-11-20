var x{1..4} integer >= -3 <= 2;

minimize o: -x[1] - 3*x[2] + x[3] - 5*x[4];

s.t. c1: (if x[1] >= 2 then x[2] else x[3]) != -2;
s.t. c2: (x[1] > -2) or (x[4] <= -3) or (x[1] >= 0);
