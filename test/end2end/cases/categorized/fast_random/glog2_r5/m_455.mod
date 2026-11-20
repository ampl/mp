var x{1..5} integer >= 0 <= 3;

maximize o: 2*x[1] - 3*x[2] + 4*x[3] - 5*x[4] + 0*x[5];

s.t. c1: (x[5] <= 3) and (x[3] >= 3) ==> (x[5] >= 3);
s.t. c2: atmost 2 {i in 1..5} (x[i] >= 1);
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
s.t. c4: ((if x[5] > 3 then 1 else 0) + (if x[5] >= 3 then 1 else 0)) >= 1;
s.t. c5: max(x[2], x[3]) - min(x[2], x[2]) <= 3;
