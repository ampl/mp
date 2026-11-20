var x{1..5} integer >= -2 <= 2;
maximize o: (-5)*x[1] + (-4)*x[2] + (4)*x[3] + (-5)*x[4] + (3)*x[5];
s.t. c1: ((if x[3] <= -2 then 1 else 0) + (if x[3] >= -2 then 1 else 0)) >= 1;
s.t. c2: (x[4] != 2) or (x[3] <= 1) or (x[4] > 0);
