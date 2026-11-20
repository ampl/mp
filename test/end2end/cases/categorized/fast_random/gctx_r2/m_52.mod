var x{1..3} integer >= 0 <= 2;
maximize o: (1)*x[1] + (-3)*x[2] + (3)*x[3];
s.t. c1: (x[1] >= x[1]) <==> (min(x[2], x[1]) <= 2);
s.t. c2: (x[1] >= x[1]) ==> ((min(x[2], x[1]) <= 2) or (abs(x[2] - x[2]) <= 1));
s.t. c3: (if abs(x[2] - x[2]) <= 1 then 1 else 0) + (if abs(x[2] - x[2]) <= 1 then 1 else 0) + (if x[1] >= x[1] then 1 else 0) >= 1;
