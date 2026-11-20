var x{1..5} integer >= 0 <= 2;
maximize o: (3)*x[1] + (4)*x[2] + (-3)*x[3] + (3)*x[4] + (2)*x[5];
s.t. c1: (min(x[3], x[3]) == 1) or (max(x[5], x[5]) != 0);
s.t. c2: not ((max(x[5], x[5]) != 0) and (x[5] == 2));
s.t. c3: (min(x[3], x[3]) == 1) ==> ((max(x[5], x[5]) != 0) or (max(x[5], x[5]) != 0));
s.t. c4: (min(x[3], x[3]) == 1) and ((max(x[5], x[5]) != 0) or not (max(x[5], x[5]) != 0));
