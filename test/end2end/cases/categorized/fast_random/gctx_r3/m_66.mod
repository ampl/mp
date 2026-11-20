var x{1..3} integer >= -1 <= 1;
maximize o: (-2)*x[1] + (0)*x[2] + (-2)*x[3];
s.t. c1: (if x[1] <= 1 then x[2] else x[3]) == 0;
s.t. c2: (x[2] + x[2] >= 1) <==> (x[2] > 0);
s.t. c3: not (x[2] > 0);
