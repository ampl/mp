var x{1..3} integer >= 0 <= 2;
maximize o: (0)*x[1] + (-4)*x[2] + (4)*x[3];
s.t. c1: abs(x[1] - x[3]) > 1;
s.t. c2: (x[1] >= 0) ==> (x[3] <= 0);
