var x{1..6} integer >= 0 <= 2;
maximize o: (6)*x[1] + (-2)*x[2] + (4)*x[3] + (5)*x[4] + (6)*x[5] + (4)*x[6];
s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 2;
s.t. c2: x[3] == 0;
s.t. c3: x[2] == 1;
