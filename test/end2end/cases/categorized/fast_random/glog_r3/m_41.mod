var x{1..3} integer >= 0 <= 2;
maximize o: (4)*x[1] + (5)*x[2] + (6)*x[3];
s.t. c1: atleast 1 {i in 1..3} (x[i] >= 2);
s.t. c2: max(x[2], x[1]) - min(x[3], x[1]) > 0;
