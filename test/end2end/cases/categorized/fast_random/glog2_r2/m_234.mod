var x{1..3} integer >= -3 <= 2;
maximize o: (5)*x[1] + (-2)*x[2] + (-1)*x[3];
s.t. c1: (x[3] < 0) or (x[3] <= -3) or (x[2] < -3);
s.t. c2: atmost 2 {i in 1..3} (x[i] >= -3);
s.t. c3: alldiff(x[2], x[1]);
