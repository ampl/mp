var x{1..3} integer >= 0 <= 2;

maximize o: -3*x[1] + 2*x[2] - 3*x[3];

s.t. c1: (x[2] <= 2) or (x[3] == 1) or (x[1] >= 2);
s.t. c2: (if x[2] <= 0 then x[1] else x[3]) >= 2;
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
