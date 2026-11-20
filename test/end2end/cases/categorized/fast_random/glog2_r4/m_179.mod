var x{1..5} integer >= 0 <= 4;

minimize o: 3*x[1] - 3*x[2] + 2*x[3] - 4*x[4] - 2*x[5];

s.t. c1: (x[4] != 3) or (x[2] == 0) or (x[4] > 0);
s.t. c2: exists {i in 1..5} x[i] == 3;
s.t. c3: exists {i in 1..5} x[i] == 1;
s.t. c4: alldiff(x[5], x[1], x[4], x[2]);
s.t. c5: (x[4] != 3) or (x[2] == 0) or (x[4] > 0);
s.t. c6: alldiff(x[5], x[1], x[4], x[2]);
