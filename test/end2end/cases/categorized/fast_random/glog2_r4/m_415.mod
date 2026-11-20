var x{1..3} integer >= 0 <= 3;
maximize o: (4)*x[1] + (6)*x[2] + (0)*x[3];
s.t. c1: (x[2] >= 2) and (x[3] == 2) ==> (x[2] >= 2);
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: exists {i in 1..3} x[i] == 3;
s.t. c4: atmost 3 {i in 1..3} (x[i] >= 3);
s.t. c5: (if x[3] == 1 then x[1] else x[2]) >= 3;
