var x{1..5} integer >= 0 <= 1;

minimize o: 6*x[1] - 5*x[2] - 4*x[3] + 3*x[4] + 2*x[5];

s.t. c1: ((if x[2] <= 0 then 1 else 0) + (if x[5] >= 0 then 1 else 0)) >= 1;
s.t. c2: (x[2] >= 1) ==> (x[5] != 0) else (x[2] != 0);
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
s.t. c4: ((if x[4] > 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
