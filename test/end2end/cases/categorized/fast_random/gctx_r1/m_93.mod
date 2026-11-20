var x{1..5} integer >= 0 <= 3;
maximize o: (0)*x[1] + (-4)*x[2] + (4)*x[3] + (-5)*x[4] + (1)*x[5] + (5)*(if min(x[1], x[5]) >= 1 then 1 else 0);
s.t. c1: (min(x[1], x[5]) >= 1) and ((min(x[1], x[1]) < 1) or not (min(x[5], x[5]) <= 0));
s.t. c2: (min(x[5], x[5]) <= 0) and ((min(x[1], x[5]) >= 1) or not (min(x[5], x[5]) <= 0));
s.t. c3: (if min(x[5], x[5]) <= 0 then x[4] else x[3]) < 3;
s.t. c4: (min(x[1], x[5]) >= 1) ==> (min(x[1], x[1]) < 1);
