var x{1..4} integer >= -2 <= 0;

minimize o: 6*x[1] - 5*x[2] + x[3] + 4*x[4];

s.t. c1: (x[2] <= -1) and (x[2] <= -2) ==> (x[2] >= -1);
s.t. c2: ((if x[1] != -1 then 1 else 0) + (if x[3] >= -1 then 1 else 0)) >= 1;
s.t. c3: (x[3] >= -1) and (x[3] >= -2) ==> (x[1] >= -1);
s.t. c4: forall {i in 1..3} x[i] != x[i+1];
