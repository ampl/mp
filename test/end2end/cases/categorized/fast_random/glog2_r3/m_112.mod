var x{1..5} integer >= 0 <= 3;
maximize o: (-3)*x[1] + (-5)*x[2] + (6)*x[3] + (0)*x[4] + (-1)*x[5];
s.t. c1: (if x[1] >= 3 then x[1] else x[1]) == 2;
s.t. c2: (if x[4] <= 3 then x[2] else x[5]) >= 3;
s.t. c3: exists {i in 1..5} x[i] == 1;
