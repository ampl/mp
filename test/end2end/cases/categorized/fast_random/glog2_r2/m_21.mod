var x{1..3} integer >= 0 <= 2;
maximize o: (1)*x[1] + (0)*x[2] + (-5)*x[3];
s.t. c1: (if x[2] <= 3 then x[1] else x[2]) <= 3;
