var x{1..5} integer >= 0 <= 3;
maximize o: (-999)*x[1] + (-3000)*x[2] + (34)*x[3] + (34)*x[4] + (1000)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: not ((x[5] == 1) and (x[1] < 2));
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
