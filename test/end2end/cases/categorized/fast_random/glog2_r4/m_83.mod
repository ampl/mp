var x{1..3} integer >= 0 <= 4;

maximize o: 6*x[1] + x[2] - x[3];

s.t. c1: max(x[3], x[3]) - min(x[3], x[3]) == 2;
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: x[3] == 2;
