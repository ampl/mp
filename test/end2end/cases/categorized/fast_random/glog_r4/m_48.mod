var x{1..4} integer >= 0 <= 3;
minimize o: (-1)*x[1] + (1)*x[2] + (4)*x[3] + (6)*x[4];
s.t. c1: ((x[2] > 1) ==> (x[1] > 3)) and ((x[1] > 3) ==> (x[1] >= 1));
s.t. c2: ((if x[4] < 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
s.t. c3: not ((x[1] > 0) and (x[4] >= 3));
