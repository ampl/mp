var x{1..5} integer >= 0 <= 4;

maximize o: -x[1] - x[2] + 6*x[3] + 6*x[4] + 0*x[5];

s.t. c1: (x[1] == 1) <==> (x[2] + x[2] < 0);
s.t. c2: not ((x[5] != 0) and (x[2] < 1));
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
s.t. c4: (x[1] >= 1) and (x[3] != 2) ==> (x[3] >= 1);
