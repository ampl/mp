var x{1..7} integer >= 0 <= 1;
maximize o: (3)*x[1] + (4)*x[2] + (4)*x[3] + (3)*x[4] + (-5)*x[5] + (5)*x[6] + (0)*x[7];
s.t. c1: (if x[4] != 1 then x[2] else x[1]) > 0;
s.t. c2: (x[4] >= 1) ==> (x[3] >= 0);
s.t. c3: abs(x[1] - x[2]) != 0;
s.t. c4: alldiff(x[4], x[7], x[6]);
