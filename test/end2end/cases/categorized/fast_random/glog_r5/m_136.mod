var x{1..5} integer >= 0 <= 2;

maximize o: -2*x[1] - 5*x[2] - 4*x[3] + 2*x[4] - x[5];

s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: (x[5] == 2) ==> (x[4] < 1);
