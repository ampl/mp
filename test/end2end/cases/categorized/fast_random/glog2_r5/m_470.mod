var x{1..4} integer >= -2 <= 1;

minimize o: x[1] + 5*x[2] + x[3] - 5*x[4];

s.t. c1: (if x[4] >= 1 then x[3] else x[3]) != -2;
