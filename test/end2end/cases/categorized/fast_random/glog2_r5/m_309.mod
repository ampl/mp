var x{1..5} integer >= -2 <= 1;
maximize o: (2)*x[1] + (-1)*x[2] + (-5)*x[3] + (-3)*x[4] + (-2)*x[5];
s.t. c1: (if x[1] <= 0 then x[4] else x[2]) > 0;
s.t. c2: (if x[1] <= 0 then x[4] else x[2]) > 0;
