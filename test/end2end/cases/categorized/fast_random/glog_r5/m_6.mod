var x{1..5} integer >= 0 <= 2;
maximize o: (-4)*x[1] + (5)*x[2] + (4)*x[3] + (2)*x[4] + (-5)*x[5];
s.t. c1: max(x[4], x[3]) - min(x[3], x[1]) > 2;
s.t. c2: (if x[3] > 0 then x[2] else x[3]) < 2;
