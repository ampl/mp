var x{1..4} integer >= 0 <= 3;
maximize o: (-4)*x[1] + (6)*x[2] + (-4)*x[3] + (-1)*x[4] + (-3)*(if x[1] == 1 then 1 else 0);
s.t. c1: (x[2] > 3) or (x[2] > 3);
s.t. c2: (x[3] != 3) or (x[2] > 3);
