var x{1..5} integer >= 0 <= 4;
maximize o: (-1)*x[1] + (5)*x[2] + (-2)*x[3] + (-2)*x[4] + (3)*x[5];
s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) == 1;
s.t. c2: (x[1] >= 2) ==> (x[4] != 1);
s.t. c3: (if x[3] >= 2 then x[5] else x[3]) >= 0;
s.t. c4: max(x[5], x[4]) - min(x[4], x[4]) < 0;
