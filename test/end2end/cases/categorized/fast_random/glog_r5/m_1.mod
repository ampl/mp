var x{1..3} integer >= -2 <= 2;

minimize o: -5*x[1] + x[2] + 6*x[3];

s.t. c1: exists {i in 1..3} x[i] == 0;
s.t. c2: (x[3] < -1) or (x[2] >= -1) or (x[1] > 0);
s.t. c3: exists {i in 1..3} x[i] == 1;
