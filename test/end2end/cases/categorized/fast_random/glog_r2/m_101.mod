var x{1..4} integer >= 0 <= 4;

minimize o: -3*x[1] - 5*x[2] + 5*x[3] + 2*x[4];

s.t. c1: exists {i in 1..4} x[i] == 3;
s.t. c2: exists {i in 1..4} x[i] == 1;
s.t. c3: abs(x[1] - x[4]) < 3;
s.t. c4: (x[1] < 4) or (x[4] < 2) or (x[4] < 0);
