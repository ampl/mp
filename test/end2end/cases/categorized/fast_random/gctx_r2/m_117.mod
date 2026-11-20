var x{1..3} integer >= 0 <= 2;
maximize o: (-1)*x[1] + (0)*x[2] + (-5)*x[3];
s.t. c1: (max(x[3], x[1]) != 2) <==> (min(x[2], x[3]) != 2);
s.t. c2: (min(x[1], x[2]) <= 1) ==> (max(x[3], x[1]) != 2) else (min(x[2], x[3]) != 2);
s.t. c3: (min(x[1], x[1]) == 1) and ((min(x[1], x[1]) == 1) or not (max(x[3], x[1]) != 2));
s.t. c4: (if min(x[1], x[2]) <= 1 then 1 else 0) + (if max(x[3], x[1]) != 2 then 1 else 0) + (if min(x[2], x[3]) != 2 then 1 else 0) >= 3;
