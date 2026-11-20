var x{1..6} integer >= 0 <= 1;
maximize o: (-5)*x[1] + (1)*x[2] + (-3)*x[3] + (5)*x[4] + (4)*x[5] + (-3)*x[6];
s.t. c1: (x[4] != 0) ==> (x[3] == 1);
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) == 1;
s.t. c3: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6]) != 2;
s.t. c4: exists {i in 1..6} x[i] == 1;
