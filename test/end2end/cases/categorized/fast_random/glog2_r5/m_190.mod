var x{1..3} integer >= 0 <= 4;
maximize o: (6)*x[1] + (2)*x[2] + (4)*x[3];
s.t. c1: numberof 3 in (x[1], x[2], x[3]) == 0;
