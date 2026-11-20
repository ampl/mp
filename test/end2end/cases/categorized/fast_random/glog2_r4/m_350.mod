var x{1..3} integer >= -2 <= 0;
maximize o: (5)*x[1] + (2)*x[2] + (-4)*x[3];
s.t. c1: exists {i in 1..3} x[i] == -1;
s.t. c2: (x[1] < -1) ==> (x[2] > -1) else (x[2] >= -1);
s.t. c3: alldiff(x[2], x[3]);
