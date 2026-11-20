var x{1..3} integer >= -2 <= 1;
maximize o: (2)*x[1] + (-1)*x[2] + (5)*x[3];
s.t. c1: (if x[2] > -2 then x[3] else x[3]) <= 2;
