var x{1..3} integer >= 0 <= 2;
maximize o: (-1)*x[1] + (6)*x[2] + (4)*x[3];
s.t. c1: alldiff(x[3], x[2]);
s.t. c2: atmost 1 {i in 1..3} (x[i] >= 1);
