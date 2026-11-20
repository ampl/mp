var x{1..6} integer >= -2 <= 1;

maximize o: 0*x[1] - 5*x[2] + 2*x[3] + x[4] + 6*x[5] - x[6];

s.t. c1: max(x[3], x[2]) - min(x[4], x[3]) != 0;
s.t. c2: abs(x[3] - x[5]) <= 1;
s.t. c3: forall {i in 1..5} x[i] != x[i+1];
s.t. c4: (x[3] > 1) and (x[5] < -2) ==> (x[3] >= 1);
