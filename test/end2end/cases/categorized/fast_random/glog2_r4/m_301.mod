var x{1..3} integer >= 0 <= 2;
maximize o: (-1)*x[1] + (2)*x[2] + (2)*x[3];
s.t. c1: (if x[1] == -1 then x[2] else x[3]) < 3;
s.t. c2: (if x[3] <= 2 then x[3] else x[2]) != 0;
s.t. c3: abs(x[3] - x[2]) >= 1;
s.t. c4: max(x[3], x[3]) - min(x[2], x[3]) <= -1;
