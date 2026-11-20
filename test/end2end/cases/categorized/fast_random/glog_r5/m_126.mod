var x{1..6} integer >= -2 <= 0;
minimize o: (5)*x[1] + (6)*x[2] + (3)*x[3] + (-3)*x[4] + (1)*x[5] + (6)*x[6];
s.t. c1: max(x[5], x[6]) - min(x[1], x[5]) >= 0;
s.t. c2: alldiff(x[1], x[6]);
s.t. c3: (if x[3] >= -1 then x[1] else x[4]) >= -2;
