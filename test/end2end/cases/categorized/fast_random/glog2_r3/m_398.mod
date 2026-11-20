var x{1..6} integer >= 0 <= 1;
maximize o: (0)*x[1] + (-1)*x[2] + (-1)*x[3] + (3)*x[4] + (0)*x[5] + (4)*x[6];
s.t. c1: max(x[6], x[5]) - min(x[4], x[2]) > 0;
s.t. c2: (x[4] <= 1) ==> (x[4] != 1) else (x[2] <= 0);
s.t. c3: exists {i in 1..6} x[i] == 0;
