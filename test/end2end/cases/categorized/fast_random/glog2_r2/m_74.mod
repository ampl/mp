var x{1..5} integer >= -3 <= 2;

maximize o: -x[1] + 0*x[2] + 5*x[3] - x[4] + 6*x[5];

s.t. c1: (x[1] > 1) and (x[1] != 2) ==> (x[1] >= 1);
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
