var x{1..4} integer >= -2 <= 3;

maximize o: 2*x[1] - 2*x[2] + 2*x[3] + 5*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: alldiff(x[4], x[3], x[2], x[1]);
