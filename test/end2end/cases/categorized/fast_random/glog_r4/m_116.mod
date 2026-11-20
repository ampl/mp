var x{1..6} integer >= -2 <= 2;

maximize o: 3*x[1] + 6*x[2] - 2*x[3] - 4*x[4] + 5*x[5] - 2*x[6];

s.t. c1: ((x[4] <= 1) ==> (x[1] > 1)) and ((x[1] > 1) ==> (x[6] >= 1));
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
