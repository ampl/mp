var x{1..3} integer >= 0 <= 4;

maximize o: 6*x[1] + 6*x[2] + x[3];

s.t. c1: abs(x[3] - x[1]) == 4;
s.t. c2: not ((x[3] == 0) and (x[3] <= 4));
s.t. c3: max(x[2], x[3]) - min(x[3], x[2]) != 2;
s.t. c4: forall {i in 1..2} x[i] != x[i+1];
