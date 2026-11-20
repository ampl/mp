var x{1..3} integer >= 0 <= 3;
maximize o: (6)*x[1] + (0)*x[2] + (-3)*x[3];
s.t. c1: (if x[1] > 1 then x[1] else x[3]) == 1;
s.t. c2: exists {i in 1..3} x[i] == 0;
