var x{1..3} integer >= 0 <= 2;

maximize o: 6*x[1] + 4*x[2] + 6*x[3];

s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: (x[2] != 1) ==> (x[3] >= 0);
s.t. c3: not ((x[2] < 2) and (x[3] == 2));
