var x{1..3} integer >= -1 <= 1;
maximize o: (3)*x[1] + (6)*x[2] + (5)*x[3];
s.t. c1: (if x[2] == 0 then x[1] else x[3]) < 0;
s.t. c2: (if x[2] == 0 then x[3] else x[3]) < 1;
