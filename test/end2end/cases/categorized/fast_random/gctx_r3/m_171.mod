var x{1..5} integer >= -1 <= 1;
maximize o: (1)*x[1] + (6)*x[2] + (-5)*x[3] + (5)*x[4] + (4)*x[5] + (3)*(if x[5] < -1 then 1 else 0);
s.t. c1: (x[1] < -1) and ((x[5] < -1) or not (x[5] < -1));
s.t. c2: (if x[1] < -1 then x[4] else 0) + (if x[1] < -1 then x[5] else 0) >= 1;
s.t. c3: (max(x[3], x[5]) > 0) ==> ((max(x[3], x[5]) > 0) or (max(x[3], x[5]) > 0));
