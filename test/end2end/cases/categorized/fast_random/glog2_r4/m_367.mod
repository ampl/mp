var x{1..6} integer >= 0 <= 3;
minimize o: (-1000)*x[1] + (2000)*x[2] + (2000)*x[3] + (2000)*x[4] + (1000)*x[5] + (2000)*x[6];
s.t. c1: (x[1] <= 1) ==> (x[3] == 0);
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: ((if x[4] == 1 then 1 else 0) + (if x[6] >= 1 then 1 else 0)) >= 1;
s.t. c4: atmost 1 {i in 1..6} (x[i] >= 1);
