var x{1..4} integer >= 0 <= 1;
maximize o: (4)*x[1] + (-4)*x[2] + (-3)*x[3] + (-5)*x[4];
s.t. c1: exists {i in 1..4} x[i] == 0;
s.t. c2: exists {i in 1..4} x[i] == 1;
s.t. c3: (x[2] != 1) ==> (x[3] > 1) else (x[2] >= 1);
s.t. c4: ((x[4] == 1) ==> (x[3] >= 1)) and ((x[3] >= 1) ==> (x[4] >= 1));
