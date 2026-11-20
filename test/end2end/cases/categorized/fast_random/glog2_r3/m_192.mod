var x{1..4} integer >= -2 <= 0;
maximize o: (3)*x[1] + (-3)*x[2] + (3)*x[3] + (5)*x[4];
s.t. c1: exists {i in 1..4} x[i] == -1;
s.t. c2: (if x[1] >= -1 then x[1] else x[1]) >= 0;
s.t. c3: (x[3] < -1) and (x[3] > 0) ==> (x[3] >= -1);
