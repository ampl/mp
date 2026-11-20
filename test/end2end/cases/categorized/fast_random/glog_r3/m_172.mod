var x{1..4} integer >= -2 <= 0;
minimize o: (4)*x[1] + (-2)*x[2] + (-3)*x[3] + (5)*x[4];
s.t. c1: not ((x[4] != -2) and (x[1] == 0));
s.t. c2: atmost 3 {i in 1..4} (x[i] >= 0);
