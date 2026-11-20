var x{1..3} integer >= 0 <= 2;

minimize o: 4*x[1] - 3*x[2] + 6*x[3];

s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: (x[1] == 2) or (x[1] <= 2) or (x[2] < 2);
s.t. c3: (x[1] > 1) and (x[2] >= 1) ==> (x[2] >= 1);
s.t. c4: x[1] == 1;
