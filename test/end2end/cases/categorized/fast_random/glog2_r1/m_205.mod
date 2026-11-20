var x{1..3} integer >= 0 <= 2;

maximize o: 3*x[1] - x[2] - x[3];

s.t. c1: exists {i in 1..3} x[i] == 0;
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: (if x[2] <= 2 then x[2] else x[1]) <= 0;
s.t. c4: abs(x[2] - x[2]) > 2;
s.t. c5: forall {i in 1..2} x[i] != x[i+1];
