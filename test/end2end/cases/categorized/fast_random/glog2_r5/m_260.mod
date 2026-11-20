var x{1..4} integer >= 0 <= 2;

maximize o: 5*x[1] + 2*x[2] + 3*x[3] - 2*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: ((if x[4] == 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
s.t. c4: numberof 0 in (x[1], x[2], x[3], x[4]) >= 3;
