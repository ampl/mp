var x{1..7} integer >= 0 <= 1;
maximize o: (6)*x[1] + (0)*x[2] + (3)*x[3] + (-1)*x[4] + (3)*x[5] + (0)*x[6] + (6)*x[7];
s.t. c1: alldiff(x[1], x[6], x[3], x[5]);
s.t. c2: (if x[3] != 0 then x[1] else x[3]) < 1;
s.t. c3: alldiff(x[2], x[6], x[1]);
