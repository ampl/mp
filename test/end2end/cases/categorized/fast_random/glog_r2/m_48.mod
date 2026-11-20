var x{1..5} integer >= -2 <= 2;
maximize o: (1)*x[1] + (5)*x[2] + (1)*x[3] + (-1)*x[4] + (5)*x[5];
s.t. c1: max(x[5], x[2]) - min(x[3], x[2]) != 1;
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: not ((x[5] >= -1) and (x[4] >= -1));
s.t. c4: (x[5] < 2) or (x[1] > -1) or (x[4] >= 0);
