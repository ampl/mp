var x{1..6} integer >= -2 <= 2;
maximize o: (-3)*x[1] + (6)*x[2] + (3)*x[3] + (4)*x[4] + (3)*x[5] + (4)*x[6];
s.t. c1: numberof -1 in (x[1], x[2], x[3], x[4], x[5], x[6]) == 3;
s.t. c2: (if x[3] > -1 then x[5] else x[6]) != 1;
s.t. c3: not ((x[3] <= -1) and (x[5] <= 2));
s.t. c4: (x[1] != 0) ==> (x[5] < 2) else (x[2] >= 2);
