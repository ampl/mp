var x{1..4} integer >= 0 <= 2;
maximize o: (2)*x[1] + (-1)*x[2] + (-5)*x[3] + (5)*x[4];
s.t. c1: (if x[1] >= 1 then x[2] else x[2]) <= 1;
