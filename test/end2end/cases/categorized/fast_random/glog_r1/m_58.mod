var x{1..4} integer >= 0 <= 4;

minimize o: -3*x[1] - 4*x[2] - 5*x[3] + 5*x[4];

s.t. c1: (x[1] >= 2) or (x[2] < 0) or (x[4] < 2);
s.t. c2: exists {i in 1..4} x[i] == 1;
