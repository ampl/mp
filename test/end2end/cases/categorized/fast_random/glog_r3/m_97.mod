var x{1..3} integer >= 0 <= 4;

maximize o: x[1] - 4*x[2] + 6*x[3];

s.t. c1: (x[1] == 1) ==> (x[2] >= 3);
s.t. c2: (x[2] == 1) ==> (x[1] > 3);
s.t. c3: not ((x[2] < 3) and (x[1] > 0));
s.t. c4: forall {i in 1..2} x[i] != x[i+1];
