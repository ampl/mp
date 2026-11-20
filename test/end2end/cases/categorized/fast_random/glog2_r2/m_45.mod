var x{1..3} integer >= -2 <= 0;
maximize o: (-5)*x[1] + (-2)*x[2] + (1)*x[3];
s.t. c1: atmost 1 {i in 1..3} (x[i] >= -1);
s.t. c2: alldiff(x[3], x[1]);
s.t. c3: (x[3] <= -1) ==> (x[3] >= -1);
s.t. c4: (x[3] >= -2) or (x[3] <= -2) or (x[3] >= 0);
s.t. c5: (x[1] <= -1) or (x[1] < 0) or (x[1] == 0);
s.t. c6: exists {i in 1..3} x[i] == -1;
s.t. c7: alldiff(x[2], x[1], x[3]);
