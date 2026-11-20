var x{1..5} integer >= -2 <= 0;
maximize o: (3)*x[1] + (-5)*x[2] + (5)*x[3] + (-1)*x[4] + (4)*x[5];
s.t. c1: not ((x[3] == -1) and (x[5] < 0));
s.t. c2: ((if x[1] == -1 then 1 else 0) + (if x[3] >= -1 then 1 else 0)) >= 1;
s.t. c3: (x[3] != 0) ==> (x[5] == 0) else (x[1] >= 0);
s.t. c4: max(x[5], x[2]) - min(x[5], x[2]) > 0;
