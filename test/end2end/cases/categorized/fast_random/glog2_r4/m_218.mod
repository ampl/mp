var x{1..5} integer >= 0 <= 1;

minimize o: -2*x[1] + 5*x[2] + 0*x[3] - 4*x[4] - 4*x[5];

s.t. c1: alldiff(x[1], x[4], x[3], x[5]);
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: (x[3] >= 1) ==> (x[3] <= 0);
s.t. c4: exists {i in 1..5} x[i] == 0;
