var x{1..5} integer >= -2 <= 2;

minimize o: 6*x[1] + 6*x[2] + x[3] - 4*x[4] - 5*x[5];

s.t. c1: (if x[1] > 0 then x[3] else x[1]) > 0;
s.t. c2: max(x[2], x[3]) - min(x[4], x[1]) > 1;
