var x{1..3} integer >= -3 <= 2;
maximize o: (2)*x[1] + (6)*x[2] + (5)*x[3];
s.t. c1: (x[3] <= -3) and (x[1] > 1) ==> (x[1] >= -3);
s.t. c2: (if x[2] < 2 then x[2] else x[3]) >= 2;
s.t. c3: exists {i in 1..3} x[i] == 1;
s.t. c4: numberof -3 in (x[1], x[2], x[3]) == 0;
