var x{1..5} integer >= -3 <= 1;
minimize o: (5)*x[1] + (0)*x[2] + (-3)*x[3] + (0)*x[4] + (5)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 0;
s.t. c2: abs(x[3] - x[4]) >= 3;
s.t. c3: max(x[4], x[3]) - min(x[5], x[2]) < 0;
