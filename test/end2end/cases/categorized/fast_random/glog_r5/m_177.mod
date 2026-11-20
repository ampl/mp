var x{1..6} integer >= 0 <= 4;
maximize o: (3)*x[1] + (2)*x[2] + (-1)*x[3] + (-3)*x[4] + (1)*x[5] + (4)*x[6];
s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 0;
