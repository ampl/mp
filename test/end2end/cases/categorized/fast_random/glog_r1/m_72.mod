var x{1..4} integer >= -2 <= 2;

maximize o: 2*x[1] - 4*x[2] - 3*x[3] - 3*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
