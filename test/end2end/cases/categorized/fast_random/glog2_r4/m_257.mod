var x{1..6} integer >= 0 <= 4;

maximize o: 2000*x[1] + 3000*x[2] + 17*x[3] - 1000*x[4] + 17*x[5] - 2000*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: atleast 2 {i in 1..6} (x[i] >= 0);
