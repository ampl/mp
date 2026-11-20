var x{1..4} integer >= -2 <= 1;

maximize o: -2*x[1] + 2*x[2] + 4*x[3] - 2*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: x[4] == 1;
s.t. c3: x[2] == -1;
