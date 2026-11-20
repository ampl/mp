var x{1..6} integer >= -2 <= 1;
minimize o: (6)*x[1] + (0)*x[2] + (-1)*x[3] + (-4)*x[4] + (3)*x[5] + (1)*x[6];
s.t. c1: atmost 2 {i in 1..6} (x[i] >= -2);
s.t. c2: not ((x[1] == -1) and (x[6] == -2));
