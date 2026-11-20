var x{1..4} integer >= -2 <= 0;
maximize o: (-5)*x[1] + (0)*x[2] + (0)*x[3] + (2)*x[4];
s.t. c1: (x[3] >= 0) and (x[1] < 0) ==> (x[4] >= 0);
s.t. c2: (if x[3] >= -2 then x[4] else x[4]) > -1;
