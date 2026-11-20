var x{1..5} integer >= 0 <= 3;

maximize o: 2*x[1] + 4*x[2] + x[3] + 5*x[4] + 6*x[5];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) == 3;
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: (x[3] >= 3) ==> (x[2] != 0) else (x[3] == 0);
s.t. c4: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) == 3;
s.t. c5: (x[3] >= 3) ==> (x[2] != 0) else (x[3] == 0);
