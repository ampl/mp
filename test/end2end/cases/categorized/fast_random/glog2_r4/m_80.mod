var x{1..6} integer >= 0 <= 3;

maximize o: 6*x[1] + 6*x[2] + 4*x[3] - 5*x[4] + 4*x[5] + 0*x[6];

s.t. c1: (x[2] <= 2) ==> (x[6] > 3);
s.t. c2: alldiff(x[3], x[4]);
s.t. c3: forall {i in 1..5} x[i] != x[i+1];
