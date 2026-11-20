var x{1..5} integer >= -1 <= 1;
maximize o: (-5)*x[1] + (-2)*x[2] + (2)*x[3] + (1)*x[4] + (0)*x[5];
s.t. c1: (if min(x[4], x[1]) != -1 then x[1] else 0) + (if min(x[4], x[1]) != -1 then x[3] else 0) < 1;
s.t. c2: (if x[1] + x[1] < 1 then x[2] else x[1]) == 0;
