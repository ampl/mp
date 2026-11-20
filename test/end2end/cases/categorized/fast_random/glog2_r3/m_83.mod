var x{1..4} integer >= 0 <= 2;
minimize o: (0)*x[1] + (-3)*x[2] + (-4)*x[3] + (1)*x[4];
s.t. c1: max(x[2], x[1]) - min(x[1], x[1]) >= 2;
s.t. c2: exists {i in 1..4} x[i] == 0;
s.t. c3: x[3] == 1;
s.t. c4: x[1] == 1;
