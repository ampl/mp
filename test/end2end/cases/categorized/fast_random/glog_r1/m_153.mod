var x{1..4} integer >= -2 <= 0;
minimize o: (-1)*x[1] + (6)*x[2] + (-3)*x[3] + (-5)*x[4];
s.t. c1: (if x[4] < -1 then x[1] else x[3]) == -1;
s.t. c2: alldiff(x[4], x[2], x[3]);
s.t. c3: (if x[4] != -1 then x[4] else x[3]) <= -1;
