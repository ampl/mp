var x{1..3} integer >= -2 <= 2;

minimize o: -x[1] + 3*x[2] + 5*x[3];

s.t. c1: (if x[3] > 1 then x[3] else x[2]) >= 0;
s.t. c2: (x[2] >= 0) ==> (x[2] != 1) else (x[1] <= 1);
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
s.t. c4: (if x[1] >= 0 then x[3] else x[1]) == 2;
