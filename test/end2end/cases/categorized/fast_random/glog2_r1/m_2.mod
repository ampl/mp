var x{1..6} integer >= -2 <= 0;
minimize o: (6)*x[1] + (-2)*x[2] + (5)*x[3] + (6)*x[4] + (6)*x[5] + (6)*x[6];
s.t. c1: (if x[3] >= -1 then x[3] else x[3]) != -2;
s.t. c2: abs(x[3] - x[3]) >= 1;
s.t. c3: (if x[1] == -1 then x[1] else x[2]) != 0;
s.t. c4: exists {i in 1..6} x[i] == -1;
