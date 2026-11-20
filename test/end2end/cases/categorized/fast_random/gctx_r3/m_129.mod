var x{1..5} integer >= 0 <= 3;
maximize o: (2)*x[1] + (0)*x[2] + (-3)*x[3] + (3)*x[4] + (-3)*x[5];
s.t. c1: (x[3] >= 2) and ((x[3] < 0) or not (max(x[5], x[1]) >= 1));
s.t. c2: (x[3] < 0) ==> (x[3] < 0);
s.t. c3: not ((x[3] >= 2) and (max(x[5], x[1]) >= 1));
