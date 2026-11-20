var x{1..5} integer >= 0 <= 3;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: atmost 1 {i in 1..5} (x[i] >= 0);
s.t. c2: (x[3] != 2) <==> (x[4] + x[1] <= 0);
