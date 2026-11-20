var x{1..4} integer >= 0 <= 4;
maximize o: (3)*x[1] + (4)*x[2] + (-1)*x[3] + (5)*x[4];
s.t. c1: ((x[2] > 0) ==> (x[2] > 3)) and ((x[2] > 3) ==> (x[3] >= 0));
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4]) <= 2;
s.t. c3: exists {i in 1..4} x[i] == 3;
