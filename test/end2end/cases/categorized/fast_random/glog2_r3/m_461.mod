var x{1..6} integer >= 0 <= 3;

minimize o: 5*x[1] + 4*x[2] + 6*x[3] - 3*x[4] - 4*x[5] + 5*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: (x[3] != 3) ==> (x[4] > 3) else (x[6] != 3);
s.t. c3: forall {i in 1..5} x[i] != x[i+1];
s.t. c4: abs(x[2] - x[3]) > 2;
