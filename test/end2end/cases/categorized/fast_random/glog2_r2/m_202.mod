var x{1..4} integer >= 0 <= 4;

maximize o: x[1] - 3*x[2] + 5*x[3] - 2*x[4];

s.t. c1: not ((x[4] == 0) and (x[4] <= 0));
s.t. c2: not ((x[2] != 4) and (x[4] == 4));
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
s.t. c4: (x[3] < 2) ==> (x[3] >= 0);
