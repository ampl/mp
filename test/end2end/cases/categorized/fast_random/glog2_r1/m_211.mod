var x{1..5} integer >= 0 <= 3;

minimize o: x[1] + 6*x[2] + 4*x[3] - x[4] - 5*x[5];

s.t. c1: abs(x[5] - x[5]) == 3;
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: (x[5] >= 0) ==> (x[5] > 3);
