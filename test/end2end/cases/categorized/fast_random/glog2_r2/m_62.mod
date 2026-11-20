var x{1..6} integer >= 0 <= 2;

minimize o: x[1] - 2*x[2] - 5*x[3] - 4*x[4] - 3*x[5] + 2*x[6];

s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 0;
s.t. c2: exists {i in 1..6} x[i] == 0;
s.t. c3: alldiff(x[5], x[6], x[2]);
