var x{1..6} integer >= 0 <= 1;

minimize o: x[1] + 4*x[2] + 0*x[3] + 3*x[4] - 5*x[5] - 2*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: alldiff(x[6], x[1], x[4]);
s.t. c3: not ((x[2] == 1) and (x[2] != 0));
