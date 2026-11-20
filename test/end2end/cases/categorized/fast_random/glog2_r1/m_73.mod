var x{1..4} integer >= -2 <= 2;

minimize o: -3*x[1] - 5*x[2] + 5*x[3] + 6*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: (x[3] != -1) ==> (x[3] != -2);
s.t. c3: x[1] == -1;
