var x{1..3} integer >= 0 <= 4;

maximize o: 3*x[1] - 5*x[2] + 4*x[3];

s.t. c1: (x[1] > 4) <==> (x[2] + x[2] < 0);
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: (x[3] >= 5) ==> (x[3] > 5);
s.t. c4: (x[3] <= 4) ==> (x[2] >= -1) else (x[1] == 5);
