var x{1..5} integer >= -2 <= 2;
minimize o: (0)*x[1] + (5)*x[2] + (3)*x[3] + (-5)*x[4] + (-4)*x[5];
s.t. c1: abs(x[1] - x[1]) <= 1;
s.t. c2: exists {i in 1..5} x[i] == 2;
s.t. c3: atmost 1 {i in 1..5} (x[i] >= 2);
s.t. c4: exists {i in 1..5} x[i] == -2;
s.t. c5: abs(x[4] - x[4]) <= 0;
s.t. c6: (x[4] == -2) or (x[2] != -1) or (x[5] == -2);
