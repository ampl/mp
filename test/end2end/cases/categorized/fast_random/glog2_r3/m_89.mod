var x{1..4} integer >= 0 <= 2;
maximize o: (4)*x[1] + (2)*x[2] + (-3)*x[3] + (0)*x[4];
s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4]) >= 2;
s.t. c2: (x[3] >= 1) and (x[3] < 2) ==> (x[1] >= 1);
s.t. c3: numberof 2 in (x[1], x[2], x[3], x[4]) >= 2;
s.t. c4: (x[3] >= 1) and (x[3] < 2) ==> (x[1] >= 1);
