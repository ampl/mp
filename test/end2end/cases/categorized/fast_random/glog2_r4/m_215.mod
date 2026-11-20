var x{1..4} integer >= -2 <= 0;
maximize o: (1)*x[1] + (-4)*x[2] + (6)*x[3] + (6)*x[4];
s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4]) >= 3;
s.t. c2: not ((x[4] != -2) and (x[1] <= -1));
s.t. c3: numberof -2 in (x[1], x[2], x[3], x[4]) <= 2;
s.t. c4: alldiff(x[4], x[3], x[2], x[1]);
s.t. c5: (x[2] <= -1) or (x[4] > -2) or (x[3] == -2);
s.t. c6: not ((x[2] == -2) and (x[2] < -1));
s.t. c7: (x[3] < 0) ==> (x[3] <= 0);
