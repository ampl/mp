var x{1..3} integer >= 0 <= 2;
maximize o: (5)*x[1] + (2)*x[2] + (-3)*x[3];
s.t. c1: (x[3] > 2) ==> (x[3] > 2);
s.t. c2: (if min(x[3], x[3]) == 2 then x[3] else x[2]) < 0;
s.t. c3: not (x[3] > 2);
s.t. c4: (x[3] > 2) ==> (min(x[3], x[3]) == 2) else (min(x[1], x[1]) <= 1);
