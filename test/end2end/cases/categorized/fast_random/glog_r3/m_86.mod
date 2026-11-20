var x{1..5} integer >= 0 <= 4;
maximize o: (-4)*x[1] + (4)*x[2] + (-3)*x[3] + (-1)*x[4] + (-2)*x[5];
s.t. c1: (if x[2] > 4 then x[4] else x[4]) == 3;
s.t. c2: (x[5] < 3) or (x[1] >= 0) or (x[3] < 3);
s.t. c3: not ((x[2] == 0) and (x[3] != 3));
s.t. c4: ((if x[5] <= 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
