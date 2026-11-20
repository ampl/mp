var x{1..4} integer >= -2 <= 1;
maximize o: (-1)*x[1] + (5)*x[2] + (-4)*x[3] + (4)*x[4];
s.t. c1: numberof -1 in (x[1], x[2], x[3], x[4]) > 0;
s.t. c2: (x[3] < 1) ==> (x[2] < -1);
