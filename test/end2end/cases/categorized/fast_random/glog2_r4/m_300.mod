var x{1..6} integer >= -2 <= 1;

maximize o: -3*x[1] - x[2] + x[3] + 4*x[4] + 6*x[5] + 2*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: not ((x[5] == 1) and (x[5] >= -2));
s.t. c3: not ((x[4] >= -1) and (x[3] < -1));
s.t. c4: not ((x[4] < 1) and (x[5] <= -1));
