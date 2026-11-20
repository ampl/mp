var x{1..6} integer >= -2 <= 0;
maximize o: (-2)*x[1] + (4)*x[2] + (3)*x[3] + (-2)*x[4] + (3)*x[5] + (0)*x[6];
s.t. c1: alldiff(x[4], x[6], x[2], x[1]);
s.t. c2: (x[4] == -1) <==> (x[4] + x[1] == -1);
s.t. c3: atmost 1 {i in 1..6} (x[i] >= -2);
