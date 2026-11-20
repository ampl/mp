var x{1..4} integer >= 0 <= 2;
minimize o: (51)*x[1] + (34)*x[2] + (-666)*x[3] + (-1000)*x[4];
s.t. c1: ((if x[4] <= 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c2: not ((x[1] < 2) and (x[3] < 0));
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
s.t. c4: (x[3] >= 0) or (x[2] < 0) or (x[3] > 1);
