var x{1..5} integer >= -2 <= 2;
maximize o: (4)*x[1] + (-2)*x[2] + (-1)*x[3] + (3)*x[4] + (6)*x[5];
s.t. c1: not ((x[2] <= -1) and (x[1] <= -2));
s.t. c2: (x[4] <= 2) ==> (x[5] < -2) else (x[5] == -1);
s.t. c3: (if x[1] >= -2 then x[3] else x[3]) != -2;
s.t. c4: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) >= 1;
s.t. c5: not ((x[2] <= -1) and (x[1] <= -2));
s.t. c6: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) >= 1;
