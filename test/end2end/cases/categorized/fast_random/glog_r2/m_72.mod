var x{1..4} integer >= 0 <= 4;

minimize o: -5*x[1] - 5*x[2] + 6*x[3] + x[4];

s.t. c1: (x[3] < 0) ==> (x[1] < 1) else (x[2] == 0);
s.t. c2: (x[1] == 2) and (x[1] >= 0) ==> (x[1] >= 2);
s.t. c3: max(x[1], x[1]) - min(x[1], x[2]) == 4;
s.t. c4: forall {i in 1..3} x[i] != x[i+1];
