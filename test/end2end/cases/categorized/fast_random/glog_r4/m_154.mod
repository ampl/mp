var x{1..6} integer >= 0 <= 4;

maximize o: -x[1] - 5*x[2] - 3*x[3] - 3*x[4] - 4*x[5] + 6*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: (x[6] == 4) ==> (x[1] >= 4) else (x[2] <= 3);
s.t. c3: abs(x[5] - x[5]) == 0;
s.t. c4: max(x[1], x[3]) - min(x[2], x[6]) != 3;
