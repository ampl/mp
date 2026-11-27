var x{1..6} integer >= 0 <= 2;

minimize o: -5*x[1] - 5*x[2] + 6*x[3] + 4*x[4] + 4*x[5] - 5*x[6];

s.t. c1: max(x[6], x[5]) - min(x[6], x[4]) != 0;
s.t. c2: abs(x[5] - x[3]) >= 2;
s.t. c3: ((x[1] != 2) ==> (x[2] > 2)) and ((x[2] > 2) ==> (x[4] >= 2));
s.t. c4: atmost 3 {i in 1..6} (x[i] >= 2);
s.t. c5: atmost 1 {i in 1..6} (x[i] >= 1);
s.t. c6: forall {i in 1..5} x[i] != x[i+1];
