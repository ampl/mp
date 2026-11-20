var x{1..4} integer >= 0 <= 3;

minimize o: 3*x[1] - 3*x[2] - 4*x[3] + 4*x[4];

s.t. c1: (if x[1] < 0 then x[3] else x[2]) <= 0;
s.t. c2: abs(x[2] - x[3]) > 2;
