var x{1..4} integer >= -2 <= 2;
minimize o: (3)*x[1] + (4)*x[2] + (0)*x[3] + (4)*x[4];
s.t. c1: alldiff(x[4], x[3], x[2]);
s.t. c2: ((if x[4] < -2 then 1 else 0) + (if x[3] >= -2 then 1 else 0)) >= 1;
s.t. c3: (x[2] >= 1) ==> (x[4] < 1) else (x[3] < -1);
s.t. c4: alldiff(x[3], x[2], x[4], x[1]);
