var x{1..6} integer >= 0 <= 2;
maximize o: (1)*x[1] + (5)*x[2] + (-5)*x[3] + (-1)*x[4] + (-2)*x[5] + (-3)*x[6];
s.t. c1: (x[3] == 0) or (x[5] != 1) or (x[6] != 1);
s.t. c2: atmost 3 {i in 1..6} (x[i] >= 0);
s.t. c3: exists {i in 1..6} x[i] == 0;
