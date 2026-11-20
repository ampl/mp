var x{1..5} integer >= 0 <= 3;

minimize o: -x[1] - x[2] + 6*x[3] + 3*x[4] - 2*x[5];

s.t. c1: (x[5] >= 1) ==> (x[1] <= 1) else (x[5] >= 0);
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: abs(x[5] - x[5]) != 2;
s.t. c4: x[1] == 0;
