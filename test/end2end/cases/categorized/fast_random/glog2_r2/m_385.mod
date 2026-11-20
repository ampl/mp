var x{1..6} integer >= 0 <= 2;
maximize o: (-4)*x[1] + (1)*x[2] + (-5)*x[3] + (4)*x[4] + (0)*x[5] + (-5)*x[6];
s.t. c1: (x[2] > 1) ==> (x[3] == 0) else (x[3] < 1);
s.t. c2: exists {i in 1..6} x[i] == 0;
s.t. c3: exists {i in 1..6} x[i] == 1;
s.t. c4: alldiff(x[3], x[1], x[4], x[5]);
s.t. c5: (x[1] >= 1) or (x[3] <= 1) or (x[3] == 0);
s.t. c6: forall {i in 1..5} x[i] != x[i+1];
