var x{1..5} integer >= 0 <= 1;
maximize o: (0)*x[1] + (-3)*x[2] + (-3)*x[3] + (5)*x[4] + (3)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 1;
s.t. c2: (if x[3] >= 1 then x[5] else x[3]) > 0;
s.t. c3: (if x[3] > 0 then x[5] else x[4]) > 1;
