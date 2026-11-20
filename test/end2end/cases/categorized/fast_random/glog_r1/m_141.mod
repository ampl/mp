var x{1..4} integer >= 0 <= 3;
maximize o: (6)*x[1] + (6)*x[2] + (5)*x[3] + (5)*x[4];
s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4]) != 3;
s.t. c2: atleast 2 {i in 1..4} (x[i] >= 2);
s.t. c3: (x[2] >= 2) ==> (x[2] == 0);
s.t. c4: abs(x[2] - x[2]) <= 0;
