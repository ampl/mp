var x{1..4} integer >= 0 <= 3;

minimize o: -666*x[1] - 666*x[2] - 1000*x[3] - 999*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: ((if x[2] == 3 then 1 else 0) + (if x[4] >= 3 then 1 else 0)) >= 1;
s.t. c3: ((if x[3] == 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
s.t. c4: (x[1] > 3) and (x[2] != 0) ==> (x[4] >= 3);
