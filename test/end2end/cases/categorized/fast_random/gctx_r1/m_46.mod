var x{1..4} integer >= 0 <= 3;
maximize o: (-2)*x[1] + (-2)*x[2] + (2)*x[3] + (-1)*x[4];
s.t. c1: (x[2] != 3) ==> ((x[3] <= x[2]) or (x[2] != 3));
s.t. c2: (max(x[4], x[2]) == 1) ==> (x[2] != 3);
s.t. c3: (if x[2] != 3 then x[1] else 0) + (if x[3] <= x[2] then x[3] else 0) < 3;
s.t. c4: (if x[2] == 1 then x[3] else x[1]) <= 3;
