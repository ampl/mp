var x{1..4} integer >= 0 <= 3;

minimize o: 0*x[1] - 5*x[2] + 5*x[3] + 3*x[4];

s.t. c1: alldiff(x[3], x[2], x[4]);
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: alldiff(x[2], x[3], x[1]);
