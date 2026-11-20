var x{1..4} integer >= -2 <= 0;

maximize o: 4*x[1] + x[2] + 3*x[3] + 2*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: max(x[4], x[1]) - min(x[1], x[2]) < -2;
