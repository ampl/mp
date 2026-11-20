var x{1..3} integer >= -2 <= 2;
maximize o: (0)*x[1] + (2)*x[2] + (0)*x[3];
s.t. c1: not ((x[2] >= -1) and (x[3] < 0));
s.t. c2: (x[2] < 1) and (x[1] <= 2) ==> (x[2] >= 1);
