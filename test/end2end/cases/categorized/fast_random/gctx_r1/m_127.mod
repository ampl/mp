var x{1..3} integer >= -1 <= 1;
maximize o: (4)*x[1] + (4)*x[2] + (-2)*x[3] + (4)*(if min(x[1], x[2]) != -1 then 1 else 0);
s.t. c1: (if x[1] > 1 then x[1] else 0) + (if x[1] > 1 then x[2] else 0) > 1;
s.t. c2: (if x[3] == -1 then x[2] else x[1]) <= 0;
s.t. c3: (x[3] == -1) ==> ((x[3] == -1) or (x[1] != -1));
s.t. c4: (x[1] != -1) or (x[3] == -1);
