var x{1..4} integer >= -3 <= 1;
minimize o: (0)*x[1] + (-5)*x[2] + (-1)*x[3] + (-3)*x[4];
s.t. c1: not ((x[2] < 0) and (x[1] < -1));
s.t. c2: (x[1] == -1) <==> (x[3] + x[3] != 0);
s.t. c3: (x[2] == -1) or (x[2] > 0) or (x[1] > 0);
s.t. c4: ((if x[1] < -2 then 1 else 0) + (if x[4] >= -2 then 1 else 0)) >= 1;
