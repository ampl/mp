var x{1..3} integer >= -2 <= 1;
maximize o: (-2)*x[1] + (5)*x[2] + (3)*x[3];
s.t. c1: atmost 2 {i in 1..3} (x[i] >= 1);
s.t. c2: (x[1] >= -2) and (x[1] >= -2) ==> (x[2] >= -2);
