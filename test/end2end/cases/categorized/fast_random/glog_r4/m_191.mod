var x{1..5} integer >= 0 <= 2;
maximize o: (-2)*x[1] + (4)*x[2] + (-3)*x[3] + (-1)*x[4] + (6)*x[5];
s.t. c1: ((x[1] <= 1) ==> (x[1] == 2)) and ((x[1] == 2) ==> (x[5] >= 1));
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: atmost 3 {i in 1..5} (x[i] >= 0);
