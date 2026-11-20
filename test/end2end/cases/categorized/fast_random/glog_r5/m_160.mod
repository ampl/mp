var x{1..6} integer >= -2 <= 1;

minimize o: x[1] - 4*x[2] + 3*x[3] + x[4] + 3*x[5] - 3*x[6];

s.t. c1: alldiff(x[6], x[2], x[3]);
s.t. c2: max(x[1], x[5]) - min(x[5], x[6]) == -2;
s.t. c3: (x[2] != 1) or (x[1] > 1) or (x[3] >= -2);
