var x{1..5} integer >= 0 <= 1;
minimize o: (3)*x[1] + (-1)*x[2] + (1)*x[3] + (2)*x[4] + (5)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: (x[2] >= 0) ==> (x[2] <= 0) else (x[3] != 1);
s.t. c3: alldiff(x[4], x[5], x[3], x[1]);
s.t. c4: max(x[2], x[1]) - min(x[1], x[3]) >= 0;
