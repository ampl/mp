var x{1..3} integer >= -2 <= 1;

maximize o: 2000*x[1] + 17*x[2] + 34*x[3];

s.t. c1: (x[1] > -2) ==> (x[1] != 0) else (x[2] != -2);
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: (x[1] <= -2) and (x[1] <= 0) ==> (x[1] >= -2);
s.t. c4: (x[2] == -2) or (x[1] != -1) or (x[3] >= 0);
