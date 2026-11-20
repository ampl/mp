var x{1..4} integer >= -2 <= 2;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: (if x[4] <= -1 then x[1] else x[4]) <= -1;
s.t. c2: not ((x[4] <= 0) and (x[4] == 2));
