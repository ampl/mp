var x{1..6} integer >= 0 <= 1;

maximize o: 2*x[1] + x[2] - 4*x[3] - 5*x[4] + 4*x[5] + 4*x[6];

s.t. c1: max(x[1], x[5]) - min(x[5], x[1]) > 1;
s.t. c2: (x[1] >= 1) ==> (x[4] > 0);
s.t. c3: forall {i in 1..5} x[i] != x[i+1];
