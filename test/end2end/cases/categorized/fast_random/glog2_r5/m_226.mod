var x{1..6} integer >= 0 <= 3;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5] + (0)*x[6];
s.t. c1: (if x[5] == 0 then x[5] else x[4]) == 1;
s.t. c2: alldiff(x[5], x[3], x[6]);
s.t. c3: (x[4] >= 1) and (x[1] >= 1) ==> (x[6] >= 1);
