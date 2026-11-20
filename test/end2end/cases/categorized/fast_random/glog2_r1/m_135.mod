var x{1..4} integer >= -2 <= 0;
maximize o: (6)*x[1] + (3)*x[2] + (-4)*x[3] + (5)*x[4];
s.t. c1: (x[1] <= -2) ==> (x[3] < 0);
s.t. c2: (x[2] == 0) ==> (x[3] < -1) else (x[4] == 0);
s.t. c3: max(x[2], x[4]) - min(x[2], x[3]) != 0;
s.t. c4: numberof -1 in (x[1], x[2], x[3], x[4]) > 1;
s.t. c5: max(x[4], x[3]) - min(x[4], x[2]) > -1;
