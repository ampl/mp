var x{1..4} integer >= 0 <= 2;
maximize o: (1)*x[1] + (-2)*x[2] + (-2)*x[3] + (-1)*x[4];
s.t. c1: (if max(x[2], x[1]) == 0 then x[4] else 0) + (if max(x[2], x[1]) == 0 then x[4] else 0) <= 2;
s.t. c2: (max(x[2], x[1]) == 0) ==> ((max(x[2], x[1]) == 0) or (x[1] <= 0));
s.t. c3: (x[3] > x[1]) ==> ((max(x[2], x[1]) == 0) or (x[1] <= 0));
s.t. c4: not (x[1] < 2);
