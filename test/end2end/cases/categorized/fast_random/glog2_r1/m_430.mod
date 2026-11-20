var x{1..6} integer >= 0 <= 2;
maximize o: (3)*x[1] + (-4)*x[2] + (-3)*x[3] + (-3)*x[4] + (5)*x[5] + (-3)*x[6];
s.t. c1: alldiff(x[4], x[6], x[5], x[2]);
s.t. c2: (if x[3] >= 0 then x[6] else x[6]) >= 0;
