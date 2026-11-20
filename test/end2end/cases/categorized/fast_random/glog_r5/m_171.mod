var x{1..3} integer >= 0 <= 3;
maximize o: (3)*x[1] + (-5)*x[2] + (6)*x[3];
s.t. c1: (if x[3] < 1 then x[1] else x[3]) <= 1;
