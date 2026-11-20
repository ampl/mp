var x{1..6} integer >= -2 <= 1;

maximize o: 0*x[1] - x[2] - 3*x[3] + 4*x[4] - 5*x[5] - 4*x[6];

s.t. c1: max(x[4], x[4]) - min(x[5], x[3]) >= -2;
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: exists {i in 1..6} x[i] == -2;
s.t. c4: forall {i in 1..5} x[i] != x[i+1];
