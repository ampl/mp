var x{1..4} integer >= 0 <= 1;

maximize o: 4*x[1] - 3*x[2] + 3*x[3] - 2*x[4];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4]) < 0;
s.t. c2: atmost 3 {i in 1..4} (x[i] >= 1);
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
