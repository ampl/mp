var x{1..5} integer >= 0 <= 2;

minimize o: -4*x[1] + 5*x[2] + x[3] + 5*x[4] - 5*x[5];

s.t. c1: (if x[5] == 0 then x[5] else x[3]) <= 0;
s.t. c2: max(x[4], x[1]) - min(x[4], x[4]) < 2;
