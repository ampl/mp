var x{1..4} integer >= -2 <= 0;
maximize o: (-5)*x[1] + (5)*x[2] + (2)*x[3] + (-5)*x[4];
s.t. c1: (x[2] == 1) ==> (x[2] > -3) else (x[4] > 0);
s.t. c2: max(x[3], x[1]) - min(x[3], x[2]) >= -3;
