var x{1..5} integer >= 0 <= 4;

minimize o: 2*x[1] + 6*x[2] + 4*x[3] - 2*x[4] - x[5];

s.t. c1: exists {i in 1..5} x[i] == 2;
s.t. c2: numberof 4 in (x[1], x[2], x[3], x[4], x[5]) < 0;
s.t. c3: exists {i in 1..5} x[i] == 2;
s.t. c4: numberof 4 in (x[1], x[2], x[3], x[4], x[5]) < 0;
