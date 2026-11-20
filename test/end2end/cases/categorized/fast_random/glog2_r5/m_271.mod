var x{1..4} integer >= 0 <= 2;

minimize o: x[1] - 5*x[2] - 2*x[3] - 4*x[4];

s.t. c1: (if x[3] <= 2 then x[1] else x[1]) <= 0;
s.t. c2: exists {i in 1..4} x[i] == -1;
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
s.t. c4: (x[4] > 0) ==> (x[1] < 0);
