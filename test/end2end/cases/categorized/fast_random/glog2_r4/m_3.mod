var x{1..6} integer >= 0 <= 4;

maximize o: -5*x[1] - 5*x[2] - x[3] - 5*x[4] + 4*x[5] - 2*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: forall {i in 1..5} x[i] != x[i+1];
s.t. c4: (x[1] >= 2) and (x[1] >= 0) ==> (x[6] >= 2);
s.t. c5: x[3] == 1;
s.t. c6: x[4] == 1;
