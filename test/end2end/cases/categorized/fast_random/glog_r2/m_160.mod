var x{1..6} integer >= 0 <= 4;
maximize o: (-1)*x[1] + (1)*x[2] + (3)*x[3] + (2)*x[4] + (-5)*x[5] + (-4)*x[6];
s.t. c1: ((if x[4] != 1 then 1 else 0) + (if x[5] >= 1 then 1 else 0)) >= 1;
