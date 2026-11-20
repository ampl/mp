var x{1..5} integer >= -2 <= 3;
maximize o: (6)*x[1] + (1)*x[2] + (4)*x[3] + (-5)*x[4] + (-5)*x[5];
s.t. c1: (x[3] == 2) ==> (x[2] < 2) else (x[5] <= 2);
s.t. c2: exists {i in 1..5} x[i] == -1;
