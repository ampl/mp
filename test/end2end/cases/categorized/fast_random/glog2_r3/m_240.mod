var x{1..3} integer >= 0 <= 3;

maximize o: 6*x[1] - 3*x[2] + 6*x[3];

s.t. c1: max(x[3], x[2]) - min(x[1], x[3]) != 1;
s.t. c2: ((x[3] <= 1) ==> (x[3] != 3)) and ((x[3] != 3) ==> (x[2] >= 1));
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
