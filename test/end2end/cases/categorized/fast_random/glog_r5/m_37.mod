var x{1..3} integer >= -2 <= 1;
maximize o: (-5)*x[1] + (5)*x[2] + (-1)*x[3];
s.t. c1: (x[1] <= -2) ==> (x[3] > -2);
s.t. c2: numberof 0 in (x[1], x[2], x[3]) > 3;
s.t. c3: abs(x[3] - x[3]) >= 2;
s.t. c4: (x[1] > 0) ==> (x[2] >= -2) else (x[1] >= -2);
