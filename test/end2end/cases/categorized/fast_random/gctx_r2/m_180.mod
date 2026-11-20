var x{1..5} integer >= 0 <= 3;
maximize o: (-1)*x[1] + (-3)*x[2] + (2)*x[3] + (6)*x[4] + (-3)*x[5];
s.t. c1: (max(x[3], x[5]) >= 3) and ((max(x[3], x[5]) >= 3) or not (x[3] <= x[3]));
s.t. c2: (if x[3] <= x[3] then 1 else 0) + (if max(x[3], x[5]) >= 3 then 1 else 0) + (if x[3] <= x[3] then 1 else 0) >= 1;
s.t. c3: (x[3] <= x[3]) ==> ((x[3] <= x[3]) or (x[3] <= x[3]));
