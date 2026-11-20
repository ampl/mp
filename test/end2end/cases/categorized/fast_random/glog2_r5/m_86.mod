var x{1..6} integer >= -2 <= 0;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5] + (0)*x[6];
s.t. c1: abs(x[6] - x[4]) > 1;
s.t. c2: (if x[6] < 0 then x[6] else x[2]) < -1;
s.t. c3: abs(x[5] - x[5]) != 1;
s.t. c4: exists {i in 1..6} x[i] == -2;
