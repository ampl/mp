var x{1..6} integer >= 0 <= 4;
maximize o: (3)*x[1] + (-5)*x[2] + (-3)*x[3] + (5)*x[4] + (2)*x[5] + (1)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 1;
s.t. c2: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 3;
