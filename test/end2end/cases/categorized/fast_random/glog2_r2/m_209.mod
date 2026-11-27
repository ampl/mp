var x{1..3} integer >= 0 <= 4;

maximize o: 6*x[1] + 2*x[2] + 0*x[3];

s.t. c1: (x[1] != 4) ==> (x[1] > 1) else (x[2] != 0);
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: numberof 3 in (x[1], x[2], x[3]) != 2;
s.t. c4: (x[1] != 4) ==> (x[1] > 1) else (x[2] != 0);
s.t. c5: numberof 3 in (x[1], x[2], x[3]) != 2;
