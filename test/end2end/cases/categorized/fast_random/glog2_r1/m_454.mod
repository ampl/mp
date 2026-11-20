var x{1..5} integer >= -2 <= 2;
maximize o: (4)*x[1] + (0)*x[2] + (1)*x[3] + (5)*x[4] + (1)*x[5];
s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) == 1;
s.t. c2: (if x[5] >= 2 then x[3] else x[3]) == -1;
