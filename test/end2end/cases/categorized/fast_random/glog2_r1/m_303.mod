var x{1..6} integer >= 0 <= 4;

minimize o: -2*x[1] - x[2] - 2*x[3] - x[4] - 4*x[5] - 5*x[6];

s.t. c1: exists {i in 1..6} x[i] == 2;
s.t. c2: x[1] == 4;
s.t. c3: x[2] == 4;
