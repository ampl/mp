var x{1..6} integer >= 0 <= 2;
maximize o: (1)*x[1] + (-4)*x[2] + (3)*x[3] + (0)*x[4] + (0)*x[5] + (5)*x[6];
s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) <= 0;
