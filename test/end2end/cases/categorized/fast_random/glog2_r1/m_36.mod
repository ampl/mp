var x{1..4} integer >= 0 <= 2;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: atmost 3 {i in 1..4} (x[i] >= 1);
s.t. c2: (x[3] != 1) ==> (x[1] < 2);
s.t. c3: alldiff(x[3], x[1], x[2]);
