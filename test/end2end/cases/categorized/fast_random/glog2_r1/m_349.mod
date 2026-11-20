var x{1..6} integer >= 0 <= 3;

maximize o: -2*x[1] - 3*x[2] + 4*x[3] - x[4] - 2*x[5] - 3*x[6];

s.t. c1: atmost 3 {i in 1..6} (x[i] >= 0);
s.t. c2: (x[6] < 3) ==> (x[1] < 2) else (x[1] == 3);
s.t. c3: max(x[5], x[3]) - min(x[1], x[1]) < 3;
s.t. c4: forall {i in 1..5} x[i] != x[i+1];
s.t. c5: atmost 3 {i in 1..6} (x[i] >= 0);
s.t. c6: forall {i in 1..5} x[i] != x[i+1];
