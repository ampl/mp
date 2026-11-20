var x{1..4} integer >= -1 <= 2;
maximize o: (1)*x[1] + (-2)*x[2] + (4)*x[3] + (-2)*x[4];
s.t. c1: not ((max(x[1], x[2]) >= 0) and (max(x[1], x[2]) >= 0));
s.t. c2: (if x[4] == 2 then 1 else 0) + (if max(x[1], x[2]) >= 0 then 1 else 0) + (if max(x[1], x[2]) >= 0 then 1 else 0) >= 2;
s.t. c3: not ((x[4] == 2) and (max(x[1], x[2]) >= 0));
s.t. c4: (if max(x[1], x[2]) >= 0 then x[1] else x[4]) >= 1;
