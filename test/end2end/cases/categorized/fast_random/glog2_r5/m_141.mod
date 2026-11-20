var x{1..5} integer >= 0 <= 4;

maximize o: x[1] + 0*x[2] - 5*x[3] - x[4] - 5*x[5];

s.t. c1: not ((x[2] > 0) and (x[1] == 0));
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
