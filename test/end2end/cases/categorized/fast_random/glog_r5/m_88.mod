var x{1..6} integer >= -2 <= 1;
minimize o: (6)*x[1] + (-2)*x[2] + (3)*x[3] + (2)*x[4] + (-4)*x[5] + (0)*x[6];
s.t. c1: alldiff(x[3], x[4], x[5], x[6]);
s.t. c2: ((if x[4] != -2 then 1 else 0) + (if x[5] >= -2 then 1 else 0)) >= 1;
s.t. c3: atleast 2 {i in 1..6} (x[i] >= -1);
