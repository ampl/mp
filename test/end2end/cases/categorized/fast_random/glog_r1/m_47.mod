var x{1..6} integer >= -2 <= 1;

maximize o: 3*x[1] + 5*x[2] + 2*x[3] - 5*x[4] + x[5] + x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: abs(x[4] - x[2]) >= 1;
s.t. c3: not ((x[6] > 0) and (x[2] >= 0));
