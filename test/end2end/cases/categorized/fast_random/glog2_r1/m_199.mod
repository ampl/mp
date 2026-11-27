var x{1..5} integer >= 0 <= 3;

maximize o: 3*x[1] + 2*x[2] + 6*x[3] + 0*x[4] + x[5];

s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: (x[1] != 0) and (x[3] < 2) ==> (x[4] >= 0);
s.t. c3: (if x[1] > 2 then x[4] else x[1]) > 1;
s.t. c4: max(x[1], x[1]) - min(x[5], x[4]) <= 1;
s.t. c5: forall {i in 1..4} x[i] != x[i+1];
s.t. c6: max(x[1], x[1]) - min(x[5], x[4]) <= 1;
