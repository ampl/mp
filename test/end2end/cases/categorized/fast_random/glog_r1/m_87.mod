var x{1..3} integer >= -2 <= 1;
minimize o: (5)*x[1] + (6)*x[2] + (4)*x[3];
s.t. c1: (x[3] <= 0) <==> (x[3] + x[2] > 1);
s.t. c2: (if x[2] >= 0 then x[1] else x[3]) > 1;
s.t. c3: alldiff(x[2], x[3], x[1]);
