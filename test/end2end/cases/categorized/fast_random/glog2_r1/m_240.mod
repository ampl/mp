var x{1..6} integer >= 0 <= 3;
maximize o: (0)*x[1] + (-1)*x[2] + (6)*x[3] + (-1)*x[4] + (6)*x[5] + (6)*x[6];
s.t. c1: not ((x[2] != 0) and (x[5] < 0));
s.t. c2: (x[4] == 0) ==> (x[3] > 2) else (x[2] > 2);
s.t. c3: ((if x[1] < 2 then 1 else 0) + (if x[1] >= 2 then 1 else 0)) >= 1;
s.t. c4: atleast 1 {i in 1..6} (x[i] >= 2);
