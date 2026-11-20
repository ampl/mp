var x{1..5} integer >= -2 <= 2;
maximize o: (5)*x[1] + (0)*x[2] + (0)*x[3] + (2)*x[4] + (1)*x[5];
s.t. c1: (if x[5] <= 2 then x[1] else x[5]) >= 0;
s.t. c2: x[2] == 0;
