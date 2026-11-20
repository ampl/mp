var x{1..3} integer >= 0 <= 3;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3];
s.t. c1: not ((x[1] == 2) and (x[2] <= 1));
s.t. c2: (x[1] > 2) and (x[3] <= 1) ==> (x[2] >= 2);
