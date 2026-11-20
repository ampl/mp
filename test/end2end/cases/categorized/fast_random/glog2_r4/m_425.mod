var x{1..4} integer >= 0 <= 4;
maximize o: (4)*x[1] + (3)*x[2] + (1)*x[3] + (-3)*x[4];
s.t. c1: exists {i in 1..4} x[i] == 3;
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: (x[4] < 2) or (x[1] <= 4) or (x[4] < 2);
s.t. c4: forall {i in 1..3} x[i] != x[i+1];
s.t. c5: (x[3] <= 0) ==> (x[2] == 1);
s.t. c6: alldiff(x[1], x[2], x[4], x[3]);
s.t. c7: exists {i in 1..4} x[i] == 3;
