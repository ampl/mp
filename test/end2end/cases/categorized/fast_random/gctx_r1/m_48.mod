var x{1..5} integer >= 0 <= 2;
maximize o: (1)*x[1] + (5)*x[2] + (3)*x[3] + (-2)*x[4] + (-3)*x[5] + (0)*(if x[2] != 1 then 1 else 0);
s.t. c1: (x[4] > x[2]) ==> ((x[2] != 1) or (x[5] == 2));
s.t. c2: (if x[2] != 1 then x[3] else x[4]) <= 2;
