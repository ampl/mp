var x{1..4} integer >= -2 <= 1;

maximize o: x[1] - 3*x[2] - 4*x[3] - 2*x[4];

s.t. c1: alldiff(x[3], x[2], x[4]);
s.t. c2: atmost 1 {i in 1..4} (x[i] >= 0);
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
s.t. c4: (x[1] > 1) ==> (x[3] != -1) else (x[2] == 0);
