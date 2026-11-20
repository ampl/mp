var x{1..5} integer >= 0 <= 3;

maximize o: 3*x[1] - x[2] + 2*x[3] - 2*x[4] + 4*x[5];

s.t. c1: (x[5] == 3) or (x[2] >= 1) or (x[5] != 0);
s.t. c2: (x[5] > 3) and (x[5] > 1) ==> (x[3] >= 3);
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
s.t. c4: x[4] == 2;
s.t. c5: x[1] == 0;
