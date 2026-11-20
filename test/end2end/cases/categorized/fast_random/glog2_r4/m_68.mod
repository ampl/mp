var x{1..4} integer >= 0 <= 1;

maximize o: -3*x[1] + x[2] + 0*x[3] + 4*x[4];

s.t. c1: (x[3] > 1) ==> (x[1] != 0) else (x[2] != 1);
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: atleast 1 {i in 1..4} (x[i] >= 1);
s.t. c4: (x[4] > 0) ==> (x[4] > 1) else (x[1] < 0);
