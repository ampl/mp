var x{1..3} integer >= 0 <= 3;
maximize o: (0)*x[1] + (0)*x[2] + (4)*x[3];
s.t. c1: (if x[1] < 3 then 1 else 0) + (if min(x[3], x[1]) != 3 then 1 else 0) + (if min(x[3], x[1]) != 3 then 1 else 0) >= 1;
s.t. c2: (min(x[3], x[1]) != 3) ==> (x[3] > 0) else (max(x[2], x[2]) >= 2);
s.t. c3: (x[3] > 0) ==> ((min(x[3], x[1]) != 3) or (min(x[3], x[1]) != 3));
