var x{1..4} integer >= 0 <= 3;

maximize o: -x[1] + x[2] + x[3] + 6*x[4];

s.t. c1: not ((x[4] > 2) and (x[4] < 3));
s.t. c2: numberof 2 in (x[1], x[2], x[3], x[4]) >= 1;
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
s.t. c4: (x[3] >= 2) ==> (x[1] == 3) else (x[2] < 3);
