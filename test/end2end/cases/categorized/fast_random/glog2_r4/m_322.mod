var x{1..4} integer >= 0 <= 4;

minimize o: -3*x[1] - 2*x[2] - 2*x[3] + x[4];

s.t. c1: (if x[1] != 1 then x[1] else x[1]) == 4;
s.t. c2: max(x[4], x[1]) - min(x[3], x[4]) >= 0;
s.t. c3: alldiff(x[2], x[3], x[1], x[4]);
