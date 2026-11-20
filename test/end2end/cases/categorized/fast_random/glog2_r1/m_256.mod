var x{1..4} integer >= -2 <= 2;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: not ((x[1] == 0) and (x[3] <= 1));
s.t. c2: ((if x[3] == 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c3: (x[3] < 2) <==> (x[3] + x[4] == 0);
s.t. c4: (x[1] >= 1) ==> (x[4] == 2);
