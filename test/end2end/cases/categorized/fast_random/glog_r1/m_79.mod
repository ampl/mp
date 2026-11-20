var x{1..3} integer >= 0 <= 3;
maximize o: (2)*x[1] + (6)*x[2] + (1)*x[3];
s.t. c1: atmost 1 {i in 1..3} (x[i] >= 3);
s.t. c2: (x[2] > 0) and (x[2] <= 2) ==> (x[1] >= 0);
s.t. c3: ((x[3] != 3) ==> (x[2] == 0)) and ((x[2] == 0) ==> (x[3] >= 3));
s.t. c4: exists {i in 1..3} x[i] == 1;
