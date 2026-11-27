var x{1..5} integer >= -2 <= 1;

maximize o: 6*x[1] + 2*x[2] + 2*x[3] - x[4] - 5*x[5];

s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: max(x[2], x[2]) - min(x[1], x[5]) <= -1;
s.t. c3: not ((x[5] >= 1) and (x[3] >= 1));
s.t. c4: atmost 1 {i in 1..5} (x[i] >= -2);
