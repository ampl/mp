var x{1..4} integer >= -2 <= 2;
maximize o: (0)*x[1] + (-5)*x[2] + (5)*x[3] + (0)*x[4];
s.t. c1: (if x[3] < 1 then x[3] else x[4]) < -1;
s.t. c2: (x[4] <= 2) and (x[3] >= 2) ==> (x[1] >= 2);
s.t. c3: exists {i in 1..4} x[i] == -1;
