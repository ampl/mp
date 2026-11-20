var x{1..3} integer >= -2 <= 1;

minimize o: -3*x[1] + 2*x[2] + 2*x[3];

s.t. c1: max(x[2], x[3]) - min(x[1], x[2]) == -1;
s.t. c2: (if x[1] < 0 then x[2] else x[2]) >= -1;
s.t. c3: alldiff(x[3], x[2]);
