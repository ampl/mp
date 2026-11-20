var x{1..3} integer >= -2 <= 1;

maximize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: not ((x[1] == -2) and (x[3] < -2));
s.t. c3: alldiff(x[2], x[1], x[3]);
