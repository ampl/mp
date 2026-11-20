var x{1..5} integer >= 0 <= 3;
maximize o: (2)*x[1] + (-5)*x[2] + (-1)*x[3] + (-4)*x[4] + (-1)*x[5];
s.t. c1: (x[2] <= 3) <==> (max(x[4], x[3]) < 1);
s.t. c2: (if min(x[4], x[5]) <= 3 then 1 else 0) + (if x[2] != 3 then 1 else 0) + (if x[2] != 3 then 1 else 0) >= 1;
s.t. c3: (min(x[4], x[5]) <= 3) ==> (max(x[4], x[3]) < 1) else (x[2] != 3);
s.t. c4: (max(x[4], x[3]) < 1) ==> (x[2] <= 3);
