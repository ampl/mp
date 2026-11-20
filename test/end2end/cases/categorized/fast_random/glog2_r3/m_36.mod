var x{1..6} integer >= 0 <= 2;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5] + (0)*x[6];
s.t. c1: alldiff(x[2], x[5]);
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) >= 0;
s.t. c4: not ((x[4] == 0) and (x[2] != 2));
