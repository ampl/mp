var x{1..4} integer >= 0 <= 4;

maximize o: -2*x[1] + 4*x[2] - x[3] - 5*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
