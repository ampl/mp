var x{1..4} integer >= 0 <= 1;
maximize o: (0)*x[1] + (4)*x[2] + (3)*x[3] + (3)*x[4];
s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4]) == 3;
