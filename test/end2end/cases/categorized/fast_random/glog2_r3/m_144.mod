var x{1..5} integer >= -2 <= 2;
minimize o: (-5)*x[1] + (2)*x[2] + (0)*x[3] + (-3)*x[4] + (0)*x[5];
s.t. c1: max(x[2], x[1]) - min(x[2], x[5]) <= -1;
s.t. c2: alldiff(x[4], x[3]);
s.t. c3: ((if x[4] >= 2 then 1 else 0) + (if x[5] >= 2 then 1 else 0)) >= 1;
