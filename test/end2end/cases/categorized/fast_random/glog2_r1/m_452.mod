var x{1..3} integer >= 0 <= 3;

minimize o: -x[1] - 3*x[2] + 5*x[3];

s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: (x[2] > 1) ==> (x[2] > 3);
s.t. c3: atmost 1 {i in 1..3} (x[i] >= 3);
s.t. c4: not ((x[2] < 3) and (x[1] > 0));
