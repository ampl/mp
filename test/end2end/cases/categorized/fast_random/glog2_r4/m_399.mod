var x{1..3} integer >= -2 <= 2;
maximize o: (6)*x[1] + (2)*x[2] + (6)*x[3];
s.t. c1: (x[3] == 2) ==> (x[2] <= -1);
s.t. c2: atmost 3 {i in 1..3} (x[i] >= -1);
s.t. c3: ((x[2] != 0) ==> (x[2] != 0)) and ((x[2] != 0) ==> (x[1] >= 0));
s.t. c4: forall {i in 1..2} x[i] != x[i+1];
s.t. c5: (x[3] == 2) ==> (x[2] <= -1);
s.t. c6: forall {i in 1..2} x[i] != x[i+1];
