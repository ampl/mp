var x{1..4} integer >= 0 <= 4;

minimize o: -333*x[1] - 333*x[2] - 333*x[3] - 1000*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: (x[1] != 3) ==> (x[2] <= 2);
s.t. c3: alldiff(x[2], x[1]);
