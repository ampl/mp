var x{1..4} integer >= 0 <= 2;

maximize o: 5*x[1] - x[2] + 6*x[3] - 4*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: numberof -1 in (x[1], x[2], x[3], x[4]) >= 3;
s.t. c3: exists {i in 1..4} x[i] == 0;
s.t. c4: alldiff(x[1], x[2]);
