var x{1..7} integer >= 0 <= 1;

minimize o: 4*x[1] - 4*x[2] + 4*x[3] - x[4] + 3*x[5] + x[6] - 2*x[7];

s.t. c1: exists {i in 1..7} x[i] == 0;
s.t. c2: alldiff(x[1], x[6], x[7]);
