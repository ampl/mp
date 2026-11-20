var x{1..5} integer >= 0 <= 3;
maximize o: (-3)*x[1] + (6)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: (if x[2] != 3 then x[3] else x[1]) >= 2;
s.t. c2: (if x[3] == 3 then 1 else 0) + (if x[2] != 3 then 1 else 0) + (if x[2] != 3 then 1 else 0) >= 2;
s.t. c3: (if x[3] == 3 then x[1] else x[5]) < 0;
