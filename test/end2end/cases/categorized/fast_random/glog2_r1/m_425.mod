var x{1..6} integer >= -2 <= 1;

minimize o: -3*x[1] - 3*x[2] - 5*x[3] + 0*x[4] - 3*x[5] + 3*x[6];

s.t. c1: atmost 3 {i in 1..6} (x[i] >= 1);
s.t. c2: atmost 1 {i in 1..6} (x[i] >= -2);
s.t. c3: exists {i in 1..6} x[i] == 1;
s.t. c4: (x[3] > -1) <==> (x[4] + x[4] <= -1);
s.t. c5: forall {i in 1..5} x[i] != x[i+1];
