var x{1..4} integer >= 0 <= 4;
minimize o: (2)*x[1] + (0)*x[2] + (4)*x[3] + (0)*x[4];
s.t. c1: abs(x[4] - x[1]) >= 2;
s.t. c2: atmost 1 {i in 1..4} (x[i] >= 2);
s.t. c3: (x[1] >= 0) or (x[2] < 3) or (x[3] != 0);
