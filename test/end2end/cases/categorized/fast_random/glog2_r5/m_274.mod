var x{1..3} integer >= -3 <= 1;
maximize o: (-2)*x[1] + (3)*x[2] + (1)*x[3];
s.t. c1: ((if x[3] >= -3 then 1 else 0) + (if x[3] >= -3 then 1 else 0)) >= 1;
s.t. c2: numberof 1 in (x[1], x[2], x[3]) != 0;
