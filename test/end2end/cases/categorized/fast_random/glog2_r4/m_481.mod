var x{1..6} integer >= 0 <= 2;
maximize o: (-2)*x[1] + (-4)*x[2] + (4)*x[3] + (4)*x[4] + (-5)*x[5] + (5)*x[6];
s.t. c1: ((if x[3] <= 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c2: (x[6] != -1) and (x[2] >= 0) ==> (x[6] >= -1);
