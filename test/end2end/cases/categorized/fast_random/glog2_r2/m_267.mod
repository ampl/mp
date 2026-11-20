var x{1..5} integer >= 0 <= 4;
maximize o: (-999)*x[1] + (-1000)*x[2] + (34)*x[3] + (17)*x[4] + (-333)*x[5];
s.t. c1: max(x[1], x[1]) - min(x[1], x[3]) < 1;
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
