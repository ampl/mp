var x{1..5} integer >= -3 <= 1;
maximize o: (-2)*x[1] + (6)*x[2] + (-2)*x[3] + (5)*x[4] + (-2)*x[5];
s.t. c1: (x[5] >= -1) or (x[3] != 1) or (x[4] > 1);
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: ((x[4] > 1) ==> (x[3] >= -2)) and ((x[3] >= -2) ==> (x[5] >= 1));
s.t. c4: exists {i in 1..5} x[i] == 0;
