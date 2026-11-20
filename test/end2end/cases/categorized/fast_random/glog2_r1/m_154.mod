var x{1..5} integer >= -3 <= 2;
maximize o: (-4)*x[1] + (-3)*x[2] + (-5)*x[3] + (-1)*x[4] + (-2)*x[5];
s.t. c1: (if x[1] != 2 then x[3] else x[3]) == -3;
