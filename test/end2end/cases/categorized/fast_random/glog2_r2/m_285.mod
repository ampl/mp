var x{1..3} integer >= 0 <= 2;
minimize o: (1)*x[1] + (3)*x[2] + (-1)*x[3];
s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: (x[3] == 0) and (x[2] > 2) ==> (x[1] >= 0);
s.t. c3: not ((x[1] < 2) and (x[2] >= 0));
s.t. c4: max(x[2], x[1]) - min(x[3], x[3]) < 1;
s.t. c5: forall {i in 1..2} x[i] != x[i+1];
s.t. c6: ((if x[2] > 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
