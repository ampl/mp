var x{1..5} integer >= 0 <= 3;
maximize o: (-1)*x[1] + (2)*x[2] + (4)*x[3] + (0)*x[4] + (-1)*x[5];
s.t. c1: (if min(x[4], x[3]) != 0 then x[1] else x[1]) == 0;
s.t. c2: (min(x[4], x[3]) != 0) ==> (abs(x[1] - x[3]) == 3);
