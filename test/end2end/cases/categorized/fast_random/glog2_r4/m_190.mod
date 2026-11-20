var x{1..5} integer >= 0 <= 2;
maximize o: (2)*x[1] + (1)*x[2] + (2)*x[3] + (2)*x[4] + (0)*x[5];
s.t. c1: atmost 2 {i in 1..5} (x[i] >= 2);
