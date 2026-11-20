var x{1..4} integer >= -2 <= 1;

minimize o: -3*x[1] - x[2] + 4*x[3] - 2*x[4];

s.t. c1: abs(x[2] - x[3]) < 1;
s.t. c2: exists {i in 1..4} x[i] == -2;
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
s.t. c4: (x[2] > 2) ==> (x[4] <= -3);
