var x{1..6} integer >= 0 <= 4;

minimize o: 4*x[1] + 0*x[2] + 0*x[3] + 4*x[4] + 2*x[5] - x[6];

s.t. c1: (x[1] != -1) ==> (x[4] < 5) else (x[4] == 4);
s.t. c2: (x[1] == 0) ==> (x[5] < 4) else (x[6] >= 5);
s.t. c3: forall {i in 1..5} x[i] != x[i+1];
