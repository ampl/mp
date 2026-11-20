var x{1..5} integer >= -1 <= 1;
maximize o: (-2)*x[1] + (-2)*x[2] + (-1)*x[3] + (1)*x[4] + (-3)*x[5];
s.t. c1: (x[5] >= x[5]) or (min(x[4], x[4]) > 1);
s.t. c2: (min(x[4], x[4]) > 1) <==> (min(x[4], x[4]) > 1);
s.t. c3: not ((x[5] >= x[5]) and (x[5] >= x[5]));
s.t. c4: (if x[5] >= x[5] then 1 else 0) + (if min(x[4], x[4]) > 1 then 1 else 0) + (if x[5] >= x[5] then 1 else 0) >= 1;
