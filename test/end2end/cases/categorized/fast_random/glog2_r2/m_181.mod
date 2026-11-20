var x{1..6} integer >= 0 <= 2;
minimize o: (1)*x[1] + (-5)*x[2] + (1)*x[3] + (2)*x[4] + (1)*x[5] + (-3)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: (x[5] >= 3) ==> (x[4] < 0) else (x[1] != 2);
s.t. c3: ((if x[1] < 2 then 1 else 0) + (if x[1] >= 2 then 1 else 0)) >= 1;
