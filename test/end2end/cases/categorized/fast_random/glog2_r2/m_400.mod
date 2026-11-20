var x{1..5} integer >= 0 <= 4;

minimize o: -4*x[1] + 2*x[2] + 2*x[3] - x[4] + 2*x[5];

s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: exists {i in 1..5} x[i] == 3;
