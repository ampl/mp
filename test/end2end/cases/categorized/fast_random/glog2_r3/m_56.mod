var x{1..5} integer >= 0 <= 4;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: (x[1] > 4) and (x[1] != 3) ==> (x[3] >= 4);
s.t. c2: numberof 3 in (x[1], x[2], x[3], x[4], x[5]) == 3;
