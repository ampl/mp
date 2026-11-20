var x{1..4} integer >= 0 <= 1;
maximize o: (0)*x[1] + (5)*x[2] + (2)*x[3] + (2)*x[4];
s.t. c1: (if x[2] < 0 then x[4] else x[3]) >= 0;
