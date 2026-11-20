var x{1..5} integer >= -2 <= 2;

maximize o: 6*x[1] + 4*x[2] - 5*x[3] - 5*x[4] + 4*x[5];

s.t. c1: (x[4] >= 2) and (x[3] != 1) ==> (x[2] >= 2);
s.t. c2: not ((x[5] != 2) and (x[2] > 1));
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
