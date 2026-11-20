var x{1..5} integer >= -2 <= 0;

minimize o: x[1] - 2*x[2] + 4*x[3] - 2*x[4] - 2*x[5];

s.t. c1: exists {i in 1..5} x[i] == -1;
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: max(x[5], x[2]) - min(x[5], x[3]) >= -1;
