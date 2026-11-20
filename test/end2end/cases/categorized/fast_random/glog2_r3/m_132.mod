var x{1..4} integer >= 0 <= 4;

minimize o: -5*x[1] - 3*x[2] + 4*x[3] - 2*x[4];

s.t. c1: not ((x[1] <= 3) and (x[2] != 2));
s.t. c2: exists {i in 1..4} x[i] == 0;
s.t. c3: (if x[4] == 1 then x[4] else x[4]) < 0;
s.t. c4: forall {i in 1..3} x[i] != x[i+1];
