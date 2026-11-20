var x{1..4} integer >= 0 <= 2;
maximize o: (0)*x[1] + (2)*x[2] + (-5)*x[3] + (2)*x[4];
s.t. c1: atmost 2 {i in 1..4} (x[i] >= 0);
s.t. c2: exists {i in 1..4} x[i] == 0;
s.t. c3: (x[1] > 1) <==> (x[1] + x[1] != 1);
s.t. c4: not ((x[2] == 1) and (x[2] <= 0));
