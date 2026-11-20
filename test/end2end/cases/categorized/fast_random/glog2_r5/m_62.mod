var x{1..4} integer >= 0 <= 2;
maximize o: (0)*x[1] + (3)*x[2] + (-1)*x[3] + (-5)*x[4];
s.t. c1: (x[1] >= 0) ==> (x[2] < 2) else (x[2] > 1);
