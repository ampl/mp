var x{1..3} integer >= 0 <= 2;
minimize o: (6)*x[1] + (3)*x[2] + (4)*x[3];
s.t. c1: exists {i in 1..3} x[i] == 0;
