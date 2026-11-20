var x{1..5} integer >= -2 <= 0;

maximize o: -4*x[1] + 3*x[2] - 3*x[3] - 2*x[4] - 4*x[5];

s.t. c1: exists {i in 1..5} x[i] == -1;
s.t. c2: max(x[4], x[3]) - min(x[3], x[1]) > 0;
s.t. c3: (x[5] >= -2) and (x[1] < -2) ==> (x[1] >= -2);
s.t. c4: forall {i in 1..4} x[i] != x[i+1];
s.t. c5: x[4] == -1;
