var x{1..5} integer >= 0 <= 2;

maximize o: -333*x[1] + 2000*x[2] + 1000*x[3] - 1000*x[4] + 17*x[5];

s.t. c1: (x[5] != 2) ==> (x[5] > 0);
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
