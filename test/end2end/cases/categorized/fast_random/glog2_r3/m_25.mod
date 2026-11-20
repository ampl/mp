var x{1..5} integer >= 0 <= 4;
maximize o: (-5)*x[1] + (-1)*x[2] + (-5)*x[3] + (3)*x[4] + (-4)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 1;
s.t. c2: (x[5] == 0) ==> (x[3] < 4) else (x[5] < 3);
s.t. c3: (if x[4] < 2 then x[5] else x[1]) == 1;
s.t. c4: (x[4] >= 1) ==> (x[1] >= 2) else (x[4] < 0);
s.t. c5: exists {i in 1..5} x[i] == 1;
