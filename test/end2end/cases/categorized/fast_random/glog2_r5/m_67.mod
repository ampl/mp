var x{1..6} integer >= 0 <= 4;
maximize o: (-999)*x[1] + (-666)*x[2] + (34)*x[3] + (51)*x[4] + (34)*x[5] + (2000)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: alldiff(x[2], x[6], x[4], x[5]);
