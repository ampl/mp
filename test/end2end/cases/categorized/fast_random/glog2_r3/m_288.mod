var x{1..4} integer >= 0 <= 1;
maximize o: (1)*x[1] + (4)*x[2] + (-5)*x[3] + (0)*x[4];
s.t. c1: (if x[1] != 1 then x[3] else x[4]) > 1;
s.t. c2: ((if x[4] < 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
