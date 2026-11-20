var x{1..4} integer >= 0 <= 3;
maximize o: (51)*x[1] + (-1000)*x[2] + (-2000)*x[3] + (-2000)*x[4];
s.t. c1: exists {i in 1..4} x[i] == 3;
s.t. c2: (x[1] <= 1) ==> (x[1] == 2) else (x[2] == 2);
