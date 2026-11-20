var x{1..3} integer >= 0 <= 3;

maximize o: -4*x[1] + 3*x[2] - x[3];

s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: abs(x[3] - x[2]) > 3;
s.t. c3: (x[3] <= 3) and (x[1] > 1) ==> (x[1] >= 3);
s.t. c4: x[3] == 2;
