var x{1..4} integer >= 0 <= 4;

maximize o: -333*x[1] - 2000*x[2] + 17*x[3] - 333*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: (x[4] >= 1) ==> (x[1] > 0) else (x[1] == 2);
