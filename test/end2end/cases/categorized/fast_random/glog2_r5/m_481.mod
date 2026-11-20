var x{1..3} integer >= 0 <= 3;
maximize o: (-2)*x[1] + (0)*x[2] + (2)*x[3];
s.t. c1: numberof 4 in (x[1], x[2], x[3]) == 1;
s.t. c2: numberof 3 in (x[1], x[2], x[3]) == 0;
s.t. c3: (x[3] == 0) or (x[2] <= 4) or (x[3] < 0);
s.t. c4: (x[1] != 0) ==> (x[3] <= 4);
