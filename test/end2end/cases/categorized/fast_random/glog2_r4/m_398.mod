var x{1..4} integer >= 0 <= 1;
maximize o: (-1)*x[1] + (6)*x[2] + (5)*x[3] + (-1)*x[4];
s.t. c1: not ((x[2] == 0) and (x[1] < 0));
s.t. c2: (x[4] != 0) or (x[3] > 0) or (x[3] <= 1);
s.t. c3: exists {i in 1..4} x[i] == 0;
