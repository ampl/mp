var x{1..7} integer >= 0 <= 1;
maximize o: (2)*x[1] + (-5)*x[2] + (6)*x[3] + (1)*x[4] + (-4)*x[5] + (6)*x[6] + (3)*x[7];
s.t. c1: alldiff(x[2], x[4], x[3], x[6]);
s.t. c2: (if x[1] > 0 then x[4] else x[1]) >= 1;
s.t. c3: (x[3] == 0) <==> (x[5] + x[3] <= 0);
