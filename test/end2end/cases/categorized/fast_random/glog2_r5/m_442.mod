var x{1..4} integer >= 0 <= 2;
maximize o: (5)*x[1] + (-4)*x[2] + (-5)*x[3] + (-4)*x[4];
s.t. c1: (if x[4] > 2 then x[1] else x[1]) <= 2;
