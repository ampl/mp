var x{1..4} integer >= 0 <= 2;
minimize o: (5)*x[1] + (1)*x[2] + (-1)*x[3] + (4)*x[4];
s.t. c1: abs(x[1] - x[2]) == 1;
s.t. c2: atmost 1 {i in 1..4} (x[i] >= 0);
