var x{1..4} integer >= -2 <= 1;
maximize o: (0)*x[1] + (3)*x[2] + (-3)*x[3] + (-3)*x[4];
s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4]) >= 0;
s.t. c2: (x[2] == -1) and (x[2] == 0) ==> (x[1] >= -1);
s.t. c3: exists {i in 1..4} x[i] == 1;
s.t. c4: numberof 1 in (x[1], x[2], x[3], x[4]) >= 0;
s.t. c5: exists {i in 1..4} x[i] == 1;
