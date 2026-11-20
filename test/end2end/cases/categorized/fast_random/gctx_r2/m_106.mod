var x{1..3} integer >= 0 <= 2;
minimize o: (-2)*x[1] + (0)*x[2] + (1)*x[3] + (2)*(if x[3] + x[1] < 1 then 1 else 0);
s.t. c1: not ((x[1] >= 1) and (x[1] == 2));
s.t. c2: (x[1] == 2) ==> (x[1] >= 1) else (x[3] + x[1] < 1);
s.t. c3: (x[1] >= 1) <==> (x[3] + x[1] < 1);
s.t. c4: not ((x[1] >= 1) and (x[3] + x[1] < 1));
