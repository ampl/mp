var x{1..3} integer >= 0 <= 4;
maximize o: (-999)*x[1] + (3000)*x[2] + (51)*x[3];
s.t. c1: alldiff(x[3], x[1], x[2]);
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: not ((x[1] != 0) and (x[3] >= 4));
s.t. c4: abs(x[1] - x[3]) >= 0;
