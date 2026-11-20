var x{1..4} integer >= -2 <= 0;
maximize o: (1)*x[1] + (4)*x[2] + (-1)*x[3] + (3)*x[4];
s.t. c1: not ((x[4] > 0) and (x[4] >= -2));
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: ((x[3] <= -2) ==> (x[1] == 0)) and ((x[1] == 0) ==> (x[4] >= -2));
s.t. c4: (x[1] >= -1) or (x[2] <= -2) or (x[2] != -2);
s.t. c5: alldiff(x[4], x[3], x[1], x[2]);
