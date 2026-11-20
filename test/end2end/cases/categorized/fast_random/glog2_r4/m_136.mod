var x{1..6} integer >= 0 <= 4;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5] + (0)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: ((if x[2] == 4 then 1 else 0) + (if x[2] >= 4 then 1 else 0)) >= 1;
s.t. c3: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) <= 2;
s.t. c4: ((x[1] != 1) ==> (x[6] < 3)) and ((x[6] < 3) ==> (x[6] >= 1));
