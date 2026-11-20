var x{1..3} integer >= -2 <= 1;
maximize o: (-1)*x[1] + (4)*x[2] + (-1)*x[3];
s.t. c1: (x[3] == -1) or (x[2] > -2) or (x[1] != 0);
s.t. c2: numberof 1 in (x[1], x[2], x[3]) == 3;
s.t. c3: abs(x[2] - x[1]) <= 0;
