var x{1..4} integer >= -2 <= 1;

maximize o: -3*x[1] + 3*x[2] - 5*x[3] - 5*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
