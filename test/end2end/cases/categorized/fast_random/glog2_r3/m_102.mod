var x{1..4} integer >= -2 <= 1;

minimize o: -4*x[1] + x[2] + 0*x[3] + 5*x[4];

s.t. c1: (x[4] != 1) and (x[4] <= 0) ==> (x[4] >= 1);
s.t. c2: atmost 2 {i in 1..4} (x[i] >= 0);
s.t. c3: exists {i in 1..4} x[i] == 1;
s.t. c4: forall {i in 1..3} x[i] != x[i+1];
