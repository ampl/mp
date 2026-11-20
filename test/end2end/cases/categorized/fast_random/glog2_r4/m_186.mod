var x{1..5} integer >= 0 <= 3;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: not ((x[4] < 3) and (x[2] != 1));
s.t. c2: not ((x[3] >= 1) and (x[4] > 2));
s.t. c3: ((if x[3] != 0 then 1 else 0) + (if x[5] >= 0 then 1 else 0)) >= 1;
