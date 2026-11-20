var x{1..3} integer >= 0 <= 2;
maximize o: (5)*x[1] + (-5)*x[2] + (4)*x[3];
s.t. c1: (min(x[3], x[3]) <= 2) ==> (min(x[3], x[3]) <= 2);
s.t. c2: (x[1] == x[2]) ==> ((x[1] == x[2]) or (min(x[3], x[3]) <= 2));
s.t. c3: (if min(x[3], x[3]) <= 2 then 1 else 0) + (if min(x[3], x[3]) <= 2 then 1 else 0) + (if x[1] == x[2] then 1 else 0) >= 3;
s.t. c4: (min(x[3], x[3]) <= 2) or (min(x[3], x[3]) <= 2);
