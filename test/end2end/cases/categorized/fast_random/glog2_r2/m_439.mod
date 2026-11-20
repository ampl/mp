var x{1..3} integer >= 0 <= 3;

minimize o: -3*x[1] - 2*x[2] + x[3];

s.t. c1: ((x[2] >= 0) ==> (x[3] != 1)) and ((x[3] != 1) ==> (x[2] >= 0));
s.t. c2: (x[3] > 1) ==> (x[1] > 1) else (x[3] > 3);
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
s.t. c4: ((x[2] >= 0) ==> (x[3] != 1)) and ((x[3] != 1) ==> (x[2] >= 0));
s.t. c5: forall {i in 1..2} x[i] != x[i+1];
