var x{1..3} integer >= -2 <= 0;

minimize o: -3*x[1] - x[2] + 4*x[3];

s.t. c1: atmost 2 {i in 1..3} (x[i] >= 0);
s.t. c2: not ((x[2] > 0) and (x[1] == -1));
s.t. c3: abs(x[3] - x[1]) < 0;
s.t. c4: forall {i in 1..2} x[i] != x[i+1];
