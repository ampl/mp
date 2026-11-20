var x{1..5} integer >= 0 <= 2;
minimize o: (4)*x[1] + (3)*x[2] + (4)*x[3] + (-5)*x[4] + (5)*x[5] + (-4)*(if x[2] >= x[2] then 1 else 0);
s.t. c1: (min(x[5], x[5]) > 1) and ((min(x[5], x[5]) > 1) or not (min(x[5], x[5]) > 1));
s.t. c2: (if x[2] >= x[2] then x[2] else 0) + (if min(x[5], x[5]) > 1 then x[3] else 0) < 1;
s.t. c3: (min(x[5], x[5]) > 1) <==> (x[2] >= x[2]);
s.t. c4: (x[2] >= x[2]) ==> (min(x[5], x[5]) > 1);
