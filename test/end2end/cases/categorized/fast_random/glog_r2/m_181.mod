var x{1..3} integer >= 0 <= 4;
minimize o: (0)*x[1] + (5)*x[2] + (1)*x[3];
s.t. c1: not ((x[3] <= 0) and (x[3] <= 1));
s.t. c2: abs(x[2] - x[3]) == 1;
s.t. c3: exists {i in 1..3} x[i] == 2;
