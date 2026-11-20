var x{1..3} integer >= 0 <= 3;
minimize o: (1)*x[1] + (1)*x[2] + (4)*x[3];
s.t. c1: (x[1] < 3) ==> (x[3] <= 3) else (x[2] >= 3);
s.t. c2: (x[3] >= 4) ==> (x[2] > 4);
s.t. c3: ((if x[1] == 4 then 1 else 0) + (if x[2] >= 4 then 1 else 0)) >= 1;
s.t. c4: (x[3] < 0) and (x[1] < 3) ==> (x[3] >= 0);
