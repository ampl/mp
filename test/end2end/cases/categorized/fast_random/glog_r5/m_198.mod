var x{1..6} integer >= 0 <= 2;

maximize o: 3*x[1] + 6*x[2] + 0*x[3] - 5*x[4] + 6*x[5] + x[6];

s.t. c1: (x[2] >= 0) ==> (x[3] != 0) else (x[6] > 0);
s.t. c2: (x[5] <= 0) ==> (x[5] < 0);
s.t. c3: forall {i in 1..5} x[i] != x[i+1];
s.t. c4: atmost 1 {i in 1..6} (x[i] >= 2);
