var x{1..3} integer >= 0 <= 2;
maximize o: (-2)*x[1] + (-5)*x[2] + (-1)*x[3];
s.t. c1: (if x[1] != -1 then x[1] else x[1]) > 3;
s.t. c2: (if x[3] != 2 then x[2] else x[1]) != -1;
