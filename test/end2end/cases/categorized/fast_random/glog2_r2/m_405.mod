var x{1..6} integer >= 0 <= 3;

minimize o: 5*x[1] + 2*x[2] + 6*x[3] + 5*x[4] - 3*x[5] - 3*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: ((if x[4] <= 2 then 1 else 0) + (if x[6] >= 2 then 1 else 0)) >= 1;
s.t. c3: ((x[4] <= 3) ==> (x[6] > 1)) and ((x[6] > 1) ==> (x[5] >= 3));
s.t. c4: forall {i in 1..5} x[i] != x[i+1];
s.t. c5: (x[4] >= 0) ==> (x[4] == 3);
