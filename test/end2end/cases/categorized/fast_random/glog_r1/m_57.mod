var x{1..5} integer >= 0 <= 4;
maximize o: (5)*x[1] + (4)*x[2] + (3)*x[3] + (-3)*x[4] + (-5)*x[5];
s.t. c1: ((if x[1] == 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
