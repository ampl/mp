var x{1..6} integer >= -2 <= 2;

maximize o: -3*x[1] + 2*x[2] - 5*x[3] + x[4] + 3*x[5] + x[6];

s.t. c1: abs(x[3] - x[1]) >= 3;
s.t. c2: not ((x[3] <= -1) and (x[5] != 2));
s.t. c3: forall {i in 1..5} x[i] != x[i+1];
