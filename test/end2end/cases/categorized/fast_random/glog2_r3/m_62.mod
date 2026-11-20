var x{1..5} integer >= 0 <= 3;

maximize o: x[1] + 2*x[2] + 0*x[3] - 4*x[4] + 2*x[5];

s.t. c1: ((x[4] < 3) ==> (x[4] != 0)) and ((x[4] != 0) ==> (x[4] >= 3));
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
s.t. c4: alldiff(x[1], x[2]);
