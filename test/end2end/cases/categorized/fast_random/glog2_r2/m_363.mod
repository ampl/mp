var x{1..3} integer >= 0 <= 4;
maximize o: (1)*x[1] + (-1)*x[2] + (0)*x[3];
s.t. c1: (x[2] >= 4) ==> (x[2] > 3);
s.t. c2: max(x[3], x[2]) - min(x[3], x[1]) <= 0;
s.t. c3: x[2] == 0;
