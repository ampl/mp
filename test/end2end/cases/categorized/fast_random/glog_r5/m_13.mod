var x{1..3} integer >= 0 <= 2;
maximize o: (-1)*x[1] + (1)*x[2] + (5)*x[3];
s.t. c1: (x[2] < 0) or (x[3] <= 0) or (x[1] == 0);
s.t. c2: numberof 2 in (x[1], x[2], x[3]) != 3;
s.t. c3: (x[2] <= 1) or (x[3] <= 2) or (x[2] > 1);
s.t. c4: abs(x[3] - x[3]) > 2;
