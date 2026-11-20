var x{1..4} integer >= 0 <= 4;
minimize o: (6)*x[1] + (-5)*x[2] + (-4)*x[3] + (5)*x[4];
s.t. c1: alldiff(x[2], x[1], x[3], x[4]);
s.t. c2: (x[4] <= 2) and (x[4] < 2) ==> (x[4] >= 2);
s.t. c3: (if x[3] == 1 then x[3] else x[4]) != 4;
