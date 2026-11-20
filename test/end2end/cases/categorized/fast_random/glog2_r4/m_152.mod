var x{1..5} integer >= -2 <= 1;
maximize o: (2)*x[1] + (6)*x[2] + (1)*x[3] + (3)*x[4] + (-3)*x[5];
s.t. c1: (x[5] >= 1) and (x[5] <= -2) ==> (x[1] >= 1);
s.t. c2: (x[5] >= 1) ==> (x[5] < 1) else (x[5] <= 1);
s.t. c3: alldiff(x[5], x[2], x[4], x[1]);
s.t. c4: numberof -1 in (x[1], x[2], x[3], x[4], x[5]) >= 1;
