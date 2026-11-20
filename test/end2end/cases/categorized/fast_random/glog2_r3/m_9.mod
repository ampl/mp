var x{1..6} integer >= 0 <= 3;
minimize o: (0)*x[1] + (4)*x[2] + (5)*x[3] + (1)*x[4] + (3)*x[5] + (1)*x[6];
s.t. c1: abs(x[3] - x[5]) <= 1;
s.t. c2: exists {i in 1..6} x[i] == 3;
s.t. c3: exists {i in 1..6} x[i] == 0;
s.t. c4: abs(x[4] - x[1]) != 2;
s.t. c5: abs(x[3] - x[5]) <= 1;
s.t. c6: abs(x[4] - x[1]) != 2;
