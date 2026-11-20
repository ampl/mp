var x{1..5} integer >= -2 <= 2;
maximize o: (-1)*x[1] + (-2)*x[2] + (4)*x[3] + (4)*x[4] + (-3)*x[5];
s.t. c1: (if x[1] >= -1 then x[2] else x[4]) <= 0;
s.t. c2: max(x[2], x[5]) - min(x[2], x[4]) < 2;
s.t. c3: (x[4] <= -1) ==> (x[4] < 1);
