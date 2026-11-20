var x{1..5} integer >= -2 <= 2;

maximize o: 3*x[1] + 2*x[2] + 0*x[3] - 4*x[4] + 4*x[5];

s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: atleast 3 {i in 1..5} (x[i] >= 0);
