var x{1..4} integer >= 0 <= 4;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: (x[2] == 4) <==> (x[1] + x[4] <= 1);
s.t. c2: not ((x[1] != 0) and (x[2] <= 0));
s.t. c3: ((if x[2] == 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c4: not ((x[3] <= 1) and (x[3] < 4));
