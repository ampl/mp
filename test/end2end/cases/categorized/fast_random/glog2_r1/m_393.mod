var x{1..4} integer >= 0 <= 3;
maximize o: (-1)*x[1] + (6)*x[2] + (-4)*x[3] + (6)*x[4];
s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4]) == 2;
s.t. c2: abs(x[4] - x[1]) <= 0;
s.t. c3: atmost 3 {i in 1..4} (x[i] >= 1);
s.t. c4: ((if x[4] <= 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c5: x[3] == 2;
