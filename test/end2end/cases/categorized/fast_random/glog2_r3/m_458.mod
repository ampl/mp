var x{1..5} integer >= 0 <= 1;
maximize o: (-2)*x[1] + (-3)*x[2] + (5)*x[3] + (3)*x[4] + (-4)*x[5];
s.t. c1: ((if x[1] == 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) < 3;
s.t. c3: atmost 2 {i in 1..5} (x[i] >= 0);
