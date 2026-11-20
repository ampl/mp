var x{1..4} integer >= 0 <= 2;
minimize o: (-1)*x[1] + (-4)*x[2] + (-1)*x[3] + (1)*x[4];
s.t. c1: ((if x[1] < 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
s.t. c2: atleast 2 {i in 1..4} (x[i] >= 1);
s.t. c3: ((if x[4] > 2 then 1 else 0) + (if x[4] >= 2 then 1 else 0)) >= 1;
s.t. c4: (x[2] == 1) ==> (x[2] <= 1) else (x[4] < 2);
