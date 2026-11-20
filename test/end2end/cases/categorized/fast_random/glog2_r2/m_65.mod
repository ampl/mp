var x{1..3} integer >= 0 <= 4;

maximize o: 4*x[1] + 2*x[2] + 6*x[3];

s.t. c1: (x[3] <= 4) and (x[1] <= 1) ==> (x[2] >= 4);
s.t. c2: (x[1] >= 2) ==> (x[3] == 4);
s.t. c3: numberof 1 in (x[1], x[2], x[3]) > 3;
s.t. c4: (x[1] != 2) ==> (x[1] != 3);
s.t. c5: forall {i in 1..2} x[i] != x[i+1];
s.t. c6: (x[1] <= 4) and (x[1] <= 3) ==> (x[1] >= 4);
