var x{1..3} integer >= 0 <= 4;

minimize o: 4*x[1] - 2*x[2] + 2*x[3];

s.t. c1: alldiff(x[3], x[1]);
s.t. c2: numberof 4 in (x[1], x[2], x[3]) >= 1;
s.t. c3: (x[1] > 4) ==> (x[1] != 3);
s.t. c4: (x[1] != 0) ==> (x[2] <= 1) else (x[1] > 1);
s.t. c5: (x[2] < 2) ==> (x[3] < 4);
s.t. c6: (x[2] < 4) or (x[1] < 4) or (x[2] <= 0);
s.t. c7: forall {i in 1..2} x[i] != x[i+1];
