var x{1..3} integer >= 0 <= 3;
minimize o: (3)*x[1] + (2)*x[2] + (3)*x[3];
s.t. c1: exists {i in 1..3} x[i] == 3;
s.t. c2: exists {i in 1..3} x[i] == 3;
