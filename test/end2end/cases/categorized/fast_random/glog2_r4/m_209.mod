var x{1..6} integer >= 0 <= 2;
maximize o: (-3)*x[1] + (2)*x[2] + (-2)*x[3] + (1)*x[4] + (-1)*x[5] + (3)*x[6];
s.t. c1: alldiff(x[6], x[5], x[2]);
s.t. c2: (if x[1] >= 1 then x[1] else x[2]) >= 1;
s.t. c3: (if x[6] != 0 then x[5] else x[3]) == 0;
s.t. c4: forall {i in 1..5} x[i] != x[i+1];
s.t. c5: alldiff(x[6], x[5], x[2]);
s.t. c6: forall {i in 1..5} x[i] != x[i+1];
