var x{1..3} integer >= 0 <= 4;
maximize o: (5)*x[1] + (1)*x[2] + (1)*x[3];
s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: max(x[1], x[3]) - min(x[1], x[3]) != 4;
s.t. c3: ((x[3] <= 4) ==> (x[1] == -1)) and ((x[1] == -1) ==> (x[1] >= 4));
s.t. c4: atmost 3 {i in 1..3} (x[i] >= 5);
