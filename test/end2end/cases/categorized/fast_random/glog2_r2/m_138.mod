var x{1..4} integer >= 0 <= 1;
maximize o: (-4)*x[1] + (1)*x[2] + (3)*x[3] + (-1)*x[4];
s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4]) <= 2;
s.t. c2: exists {i in 1..4} x[i] == 1;
s.t. c3: (x[3] != 0) <==> (x[1] + x[4] <= 1);
