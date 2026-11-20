var x{1..4} integer >= 0 <= 3;
maximize o: (4)*x[1] + (1)*x[2] + (2)*x[3] + (4)*x[4] + (4)*(if max(x[3], x[1]) <= 3 then 1 else 0);
s.t. c1: (if x[4] >= x[3] then x[4] else 0) + (if x[1] + x[1] != 1 then x[4] else 0) > 0;
s.t. c2: (if max(x[3], x[1]) <= 3 then 1 else 0) + (if max(x[3], x[1]) <= 3 then 1 else 0) + (if x[1] + x[1] != 1 then 1 else 0) >= 2;
s.t. c3: (x[1] + x[1] != 1) or (x[1] + x[1] != 1);
s.t. c4: (x[4] >= x[3]) ==> ((x[4] >= x[3]) or (x[4] >= x[3]));
