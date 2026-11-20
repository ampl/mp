var x{1..4} integer >= -2 <= 2;

maximize o: 3*x[1] + 0*x[2] + 5*x[3] + 5*x[4];

s.t. c1: abs(x[4] - x[2]) >= 0;
s.t. c2: not ((x[3] > 1) and (x[1] <= 0));
s.t. c3: exists {i in 1..4} x[i] == 1;
s.t. c4: forall {i in 1..3} x[i] != x[i+1];
s.t. c5: x[4] == 1;
s.t. c6: x[3] == 0;
