var x{1..3} integer >= 0 <= 2;

minimize o: 6*x[1] - 5*x[2] - 4*x[3];

s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: (x[2] > 0) or (x[3] == 0) or (x[3] <= 0);
s.t. c4: (x[3] > 0) and (x[2] >= 2) ==> (x[3] >= 0);
