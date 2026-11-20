var x{1..5} integer >= -2 <= 1;

minimize o: x[1] + 4*x[2] + x[3] - 2*x[4] - 3*x[5];

s.t. c1: alldiff(x[3], x[4], x[2]);
s.t. c2: abs(x[2] - x[5]) != 3;
s.t. c3: exists {i in 1..5} x[i] == -1;
