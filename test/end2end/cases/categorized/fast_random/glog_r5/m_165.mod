var x{1..5} integer >= 0 <= 3;
maximize o: (2)*x[1] + (4)*x[2] + (3)*x[3] + (3)*x[4] + (5)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 1;
s.t. c2: (x[3] >= 0) ==> (x[2] != 3) else (x[4] < 0);
s.t. c3: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) < 1;
