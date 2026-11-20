var x{1..4} integer >= 0 <= 2;

maximize o: -x[1] + 5*x[2] - 2*x[3] + 3*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: not ((x[1] < 2) and (x[2] != -1));
