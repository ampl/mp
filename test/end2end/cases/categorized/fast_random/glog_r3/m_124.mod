var x{1..6} integer >= 0 <= 2;

minimize o: -2*x[1] - 5*x[2] + 3*x[3] + x[4] - 2*x[5] + 6*x[6];

s.t. c1: not ((x[5] <= 0) and (x[5] >= 2));
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: (x[1] == 0) or (x[6] >= 0) or (x[6] >= 2);
s.t. c4: (x[2] == 0) ==> (x[3] >= 1);
