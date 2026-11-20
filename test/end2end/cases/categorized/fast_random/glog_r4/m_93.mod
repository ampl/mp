var x{1..3} integer >= 0 <= 3;
maximize o: (1)*x[1] + (-1)*x[2] + (-3)*x[3];
s.t. c1: numberof 3 in (x[1], x[2], x[3]) < 1;
s.t. c2: (if x[1] < 1 then x[2] else x[3]) <= 2;
s.t. c3: not ((x[3] < 2) and (x[1] < 3));
