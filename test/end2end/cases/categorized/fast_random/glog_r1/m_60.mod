var x{1..3} integer >= 0 <= 2;

minimize o: x[1] - 3*x[2] + 3*x[3];

s.t. c1: alldiff(x[2], x[1]);
s.t. c2: (x[2] < 0) or (x[2] > 2) or (x[1] != 0);
