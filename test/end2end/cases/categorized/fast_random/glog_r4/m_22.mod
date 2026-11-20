var x{1..5} integer >= 0 <= 2;
maximize o: (6)*x[1] + (5)*x[2] + (4)*x[3] + (-2)*x[4] + (3)*x[5];
s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) == 0;
s.t. c2: (x[2] < 2) ==> (x[3] > 2) else (x[3] != 2);
s.t. c3: alldiff(x[3], x[5], x[1], x[2]);
s.t. c4: (if x[5] < 2 then x[4] else x[2]) <= 0;
