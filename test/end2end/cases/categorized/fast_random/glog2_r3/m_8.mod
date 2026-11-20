var x{1..4} integer >= 0 <= 1;
maximize o: (-4)*x[1] + (4)*x[2] + (6)*x[3] + (-4)*x[4];
s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4]) <= 2;
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4]) >= 0;
s.t. c3: atmost 3 {i in 1..4} (x[i] >= 0);
