var x{1..5} integer >= 0 <= 2;
maximize o: (-2)*x[1] + (-2)*x[2] + (4)*x[3] + (2)*x[4] + (2)*x[5];
s.t. c1: (if x[2] + x[2] <= 1 then x[5] else x[2]) >= 0;
s.t. c2: (min(x[2], x[2]) > 1) ==> (min(x[2], x[2]) > 1);
s.t. c3: (min(x[2], x[2]) > 1) <==> (x[2] + x[2] <= 1);
