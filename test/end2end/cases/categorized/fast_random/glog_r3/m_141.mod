var x{1..4} integer >= -2 <= 1;
maximize o: (3)*x[1] + (5)*x[2] + (-3)*x[3] + (-5)*x[4];
s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4]) < 2;
