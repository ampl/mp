var x{1..4} integer >= 0 <= 2;

maximize o: x[1] - 2*x[2] + 3*x[3] - x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: (x[4] > 2) <==> (x[1] + x[2] < 0);
s.t. c3: max(x[2], x[2]) - min(x[3], x[2]) <= 1;
