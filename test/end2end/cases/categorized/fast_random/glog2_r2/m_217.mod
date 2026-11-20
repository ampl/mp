var x{1..4} integer >= 0 <= 3;

minimize o: -666*x[1] + 2000*x[2] - 2000*x[3] - 1000*x[4];

s.t. c1: exists {i in 1..4} x[i] == 1;
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
