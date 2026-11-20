var x{1..5} integer >= 0 <= 1;
maximize o: (-5)*x[1] + (0)*x[2] + (2)*x[3] + (0)*x[4] + (5)*x[5];
s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) <= 0;
s.t. c2: (x[4] < 1) and (x[1] >= 0) ==> (x[4] >= 1);
