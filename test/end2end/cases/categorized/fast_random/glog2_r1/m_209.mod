var x{1..6} integer >= -2 <= 2;
maximize o: (-3)*x[1] + (-2)*x[2] + (-3)*x[3] + (0)*x[4] + (-1)*x[5] + (0)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 1;
s.t. c2: (x[5] <= 0) and (x[3] > 1) ==> (x[4] >= 0);
s.t. c3: exists {i in 1..6} x[i] == 1;
s.t. c4: (x[5] <= 0) and (x[3] > 1) ==> (x[4] >= 0);
