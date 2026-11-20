var x{1..4} integer >= -3 <= 1;

minimize o: 3*x[1] - x[2] + 0*x[3] - 5*x[4];

s.t. c1: ((if x[3] < 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c2: ((x[1] > 0) ==> (x[1] <= 1)) and ((x[1] <= 1) ==> (x[1] >= 0));
s.t. c3: max(x[4], x[1]) - min(x[2], x[3]) < -1;
s.t. c4: forall {i in 1..3} x[i] != x[i+1];
