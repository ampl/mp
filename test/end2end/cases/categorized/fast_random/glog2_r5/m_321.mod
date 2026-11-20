var x{1..5} integer >= 0 <= 4;
minimize o: (-4)*x[1] + (-4)*x[2] + (1)*x[3] + (3)*x[4] + (0)*x[5];
s.t. c1: (x[4] == 5) ==> (x[4] <= 0) else (x[3] >= 4);
s.t. c2: ((if x[2] > 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c3: (if x[3] > 0 then x[2] else x[4]) == -1;
s.t. c4: atleast 3 {i in 1..5} (x[i] >= 0);
