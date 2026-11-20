var x{1..5} integer >= -1 <= 2;
maximize o: (-2)*x[1] + (0)*x[2] + (4)*x[3] + (5)*x[4] + (1)*x[5];
s.t. c1: not (min(x[2], x[2]) == 0);
s.t. c2: (x[3] <= x[1]) ==> (x[3] <= x[1]) else (x[3] <= x[1]);
