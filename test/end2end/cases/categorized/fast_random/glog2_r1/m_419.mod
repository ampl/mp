var x{1..5} integer >= -2 <= 1;
maximize o: (-3)*x[1] + (0)*x[2] + (3)*x[3] + (4)*x[4] + (0)*x[5];
s.t. c1: (if x[3] != -2 then x[1] else x[3]) == -2;
s.t. c2: not ((x[1] == -1) and (x[2] > 1));
s.t. c3: (x[5] >= -1) ==> (x[1] >= -2) else (x[4] > 1);
s.t. c4: (if x[3] != -2 then x[1] else x[3]) == -2;
s.t. c5: (x[5] >= -1) ==> (x[1] >= -2) else (x[4] > 1);
