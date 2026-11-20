var x{1..5} integer >= 0 <= 2;
maximize o: (-3)*x[1] + (5)*x[2] + (6)*x[3] + (-2)*x[4] + (-2)*x[5];
s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) >= 3;
