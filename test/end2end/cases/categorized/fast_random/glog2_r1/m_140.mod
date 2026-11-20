var x{1..3} integer >= 0 <= 4;

minimize o: 2*x[1] - 2*x[2] - 2*x[3];

s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: (x[2] == 2) <==> (x[2] + x[3] <= 4);
