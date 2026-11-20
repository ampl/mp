var x{1..4} integer >= 0 <= 2;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: (x[2] < 1) ==> (x[2] > 0) else (x[3] < 0);
