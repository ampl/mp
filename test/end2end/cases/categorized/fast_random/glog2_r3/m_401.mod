var x{1..3} integer >= 0 <= 2;
minimize o: (-3)*x[1] + (3)*x[2] + (5)*x[3];
s.t. c1: not ((x[3] == -1) and (x[1] == -1));
s.t. c2: ((if x[2] < 3 then 1 else 0) + (if x[3] >= 3 then 1 else 0)) >= 1;
s.t. c3: alldiff(x[3], x[1]);
