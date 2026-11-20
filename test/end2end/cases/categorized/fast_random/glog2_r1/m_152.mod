var x{1..5} integer >= 0 <= 3;
minimize o: (-3)*x[1] + (1)*x[2] + (3)*x[3] + (-5)*x[4] + (3)*x[5];
s.t. c1: (x[5] > 1) or (x[5] >= 3) or (x[5] <= 3);
s.t. c2: exists {i in 1..5} x[i] == 2;
s.t. c3: atmost 1 {i in 1..5} (x[i] >= 3);
