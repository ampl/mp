var x{1..4} integer >= -2 <= 0;

minimize o: x[1] + x[2] - 3*x[3] + 0*x[4];

s.t. c1: ((if x[3] < -2 then 1 else 0) + (if x[1] >= -2 then 1 else 0)) >= 1;
s.t. c2: alldiff(x[4], x[2]);
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
s.t. c4: exists {i in 1..4} x[i] == -2;
