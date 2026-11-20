var x{1..4} integer >= 0 <= 4;
minimize o: (1)*x[1] + (-4)*x[2] + (-3)*x[3] + (5)*x[4];
s.t. c1: alldiff(x[2], x[4]);
s.t. c2: (x[3] > 1) or (x[4] > 2) or (x[4] <= 2);
s.t. c3: (if x[3] > 1 then x[3] else x[1]) < 3;
