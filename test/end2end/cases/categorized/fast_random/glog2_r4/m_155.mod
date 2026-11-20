var x{1..4} integer >= -2 <= 1;

maximize o: 3*x[1] + 2*x[2] - 4*x[3] - 4*x[4];

s.t. c1: alldiff(x[2], x[4], x[1], x[3]);
s.t. c2: (x[2] != 1) or (x[4] != 0) or (x[3] != -1);
s.t. c3: numberof -1 in (x[1], x[2], x[3], x[4]) <= 2;
s.t. c4: forall {i in 1..3} x[i] != x[i+1];
s.t. c5: atmost 3 {i in 1..4} (x[i] >= 0);
